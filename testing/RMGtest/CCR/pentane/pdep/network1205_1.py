species(
    label = '[CH2][CH]C(C)C(=C)[O](2743)',
    structure = SMILES('[CH2][CH]C(C)C(=C)[O]'),
    E0 = (250.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,997.736,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0268206,'amu*angstrom^2'), symmetry=1, barrier=(18.9471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00739858,'amu*angstrom^2'), symmetry=1, barrier=(5.22675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227331,'amu*angstrom^2'), symmetry=1, barrier=(5.22678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.824074,'amu*angstrom^2'), symmetry=1, barrier=(18.9471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645949,0.0649023,-4.90165e-05,1.93488e-08,-3.09336e-12,30314.1,32.2727], Tmin=(100,'K'), Tmax=(1479.08,'K')), NASAPolynomial(coeffs=[14.9066,0.0263364,-9.90575e-06,1.7206e-09,-1.13823e-13,26095.5,-42.1091], Tmin=(1479.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(RCCJ)"""),
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
    label = 'm1_allyl(186)',
    structure = SMILES('C=C[CH]C'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800698,'amu*angstrom^2'), symmetry=1, barrier=(1.84096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.062381,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC(C)C(=C)[O](2741)',
    structure = SMILES('C=CC(C)C(=C)[O]'),
    E0 = (-21.5737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,326.269,326.427,327.172],'cm^-1')),
        HinderedRotor(inertia=(0.134119,'amu*angstrom^2'), symmetry=1, barrier=(10.1419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133961,'amu*angstrom^2'), symmetry=1, barrier=(10.1419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133922,'amu*angstrom^2'), symmetry=1, barrier=(10.1391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.565508,0.0683287,-5.57221e-05,2.41878e-08,-4.25649e-12,-2465.11,26.2997], Tmin=(100,'K'), Tmax=(1353.73,'K')), NASAPolynomial(coeffs=[14.4619,0.0272675,-1.02241e-05,1.78152e-09,-1.18598e-13,-6227.5,-44.9512], Tmin=(1353.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.5737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CH2(S)(11)',
    structure = SMILES('[CH2]'),
    E0 = (418.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1381.44,2681.05,3148.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08586,-0.00146571,5.70225e-06,-3.92867e-09,8.88752e-13,50365.3,-0.343718], Tmin=(100,'K'), Tmax=(1285.36,'K')), NASAPolynomial(coeffs=[2.65113,0.0038016,-1.38109e-06,2.30893e-10,-1.47438e-14,50667.9,6.68035], Tmin=(1285.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]CC([CH2])=O(1905)',
    structure = SMILES('[CH2][CH]CC(=C)[O]'),
    E0 = (282.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,474.042,1951.29,1951.66],'cm^-1')),
        HinderedRotor(inertia=(0.0629251,'amu*angstrom^2'), symmetry=1, barrier=(9.96073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0155193,'amu*angstrom^2'), symmetry=1, barrier=(2.70317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0637649,'amu*angstrom^2'), symmetry=1, barrier=(9.9335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.1,'J/mol'), sigma=(6.27192,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.76 K, Pc=33.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67971,0.0499013,-4.01055e-05,1.84345e-08,-3.58896e-12,34077.8,26.7309], Tmin=(100,'K'), Tmax=(1197.47,'K')), NASAPolynomial(coeffs=[8.80395,0.0261036,-1.02955e-05,1.83839e-09,-1.24118e-13,32371.6,-8.92332], Tmin=(1197.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)[CH]C(=C)[O](4200)',
    structure = SMILES('[CH2]C([O])=CC([CH2])C'),
    E0 = (202.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0738434,0.0723192,-6.36448e-05,2.97756e-08,-5.44134e-12,24562,32.0109], Tmin=(100,'K'), Tmax=(1469.69,'K')), NASAPolynomial(coeffs=[16.5266,0.0216074,-5.83176e-06,8.04262e-10,-4.59666e-14,20366.7,-51.5193], Tmin=(1469.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([CH]C)C(=C)[O](2503)',
    structure = SMILES('[CH2]C([CH]C)C(=C)[O]'),
    E0 = (250.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,1962.86],'cm^-1')),
        HinderedRotor(inertia=(0.0102752,'amu*angstrom^2'), symmetry=1, barrier=(2.7886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0102751,'amu*angstrom^2'), symmetry=1, barrier=(2.78855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373544,'amu*angstrom^2'), symmetry=1, barrier=(10.1375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373544,'amu*angstrom^2'), symmetry=1, barrier=(10.1375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586634,0.0664197,-5.34853e-05,2.33575e-08,-4.13471e-12,30296.4,32.4745], Tmin=(100,'K'), Tmax=(1352.09,'K')), NASAPolynomial(coeffs=[14.0751,0.0265156,-9.2158e-06,1.52979e-09,-9.87714e-14,26648.9,-36.6682], Tmin=(1352.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH][CH2](502)',
    structure = SMILES('[CH][CH2]'),
    E0 = (557.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1433.18,1433.5],'cm^-1')),
        HinderedRotor(inertia=(0.00559429,'amu*angstrom^2'), symmetry=1, barrier=(8.15686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7298,0.00311091,1.59755e-05,-2.26123e-08,9.32732e-12,67119.1,8.34543], Tmin=(100,'K'), Tmax=(875.199,'K')), NASAPolynomial(coeffs=[4.97018,0.00511238,-6.01253e-07,2.8794e-11,-6.01916e-16,66608.3,0.848354], Tmin=(875.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C([O])[CH]C(2381)',
    structure = SMILES('[CH2]C([O])=CC'),
    E0 = (52.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,271.139],'cm^-1')),
        HinderedRotor(inertia=(0.126869,'amu*angstrom^2'), symmetry=1, barrier=(6.61654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126818,'amu*angstrom^2'), symmetry=1, barrier=(6.61712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96216,0.0395244,-2.15034e-05,-2.39371e-09,4.75124e-12,6363.74,18.5337], Tmin=(100,'K'), Tmax=(938.593,'K')), NASAPolynomial(coeffs=[10.4571,0.0163491,-5.28597e-06,8.75321e-10,-5.83608e-14,4195.25,-24.968], Tmin=(938.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
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
    label = '[CH2][CH]C(C)[C]=C(4201)',
    structure = SMILES('[CH2][CH]C(C)[C]=C'),
    E0 = (565.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180,766.423],'cm^-1')),
        HinderedRotor(inertia=(0.00990729,'amu*angstrom^2'), symmetry=1, barrier=(4.12936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179612,'amu*angstrom^2'), symmetry=1, barrier=(4.12963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179599,'amu*angstrom^2'), symmetry=1, barrier=(4.12934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0099055,'amu*angstrom^2'), symmetry=1, barrier=(4.12883,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61024,0.0505811,-2.94689e-05,8.41683e-09,-9.82653e-13,68110.1,27.7172], Tmin=(100,'K'), Tmax=(1880.26,'K')), NASAPolynomial(coeffs=[12.216,0.0280187,-1.14694e-05,2.0349e-09,-1.34107e-13,64121.8,-30.1459], Tmin=(1880.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S)"""),
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
    label = '[CH2][C]C(C)C(=C)[O](4202)',
    structure = SMILES('[CH2][C]C(C)C(=C)[O]'),
    E0 = (504.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,354.612,354.673,354.722,354.784],'cm^-1')),
        HinderedRotor(inertia=(0.23219,'amu*angstrom^2'), symmetry=1, barrier=(20.7212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23208,'amu*angstrom^2'), symmetry=1, barrier=(20.7193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00292286,'amu*angstrom^2'), symmetry=1, barrier=(20.719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00133965,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561554,0.0674925,-5.21764e-05,1.57847e-08,-2.89493e-14,60826.7,29.913], Tmin=(100,'K'), Tmax=(1010.53,'K')), NASAPolynomial(coeffs=[15.5625,0.0230322,-8.32504e-06,1.4639e-09,-1.00163e-13,57033.2,-46.3841], Tmin=(1010.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C(C)C(=C)[O](4203)',
    structure = SMILES('[CH][CH]C(C)C(=C)[O]'),
    E0 = (493.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626373,0.0664894,-5.59305e-05,2.45652e-08,-4.33096e-12,59536.2,31.5544], Tmin=(100,'K'), Tmax=(1359.56,'K')), NASAPolynomial(coeffs=[15.0606,0.0240223,-9.0769e-06,1.59042e-09,-1.06316e-13,55611.4,-42.5162], Tmin=(1359.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1OC(=C)C1C(4204)',
    structure = SMILES('[CH2]C1OC(=C)C1C'),
    E0 = (87.0973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719462,0.0522478,2.23993e-05,-8.21151e-08,4.0703e-11,10612.2,19.6503], Tmin=(100,'K'), Tmax=(903.939,'K')), NASAPolynomial(coeffs=[21.723,0.0125936,-2.24734e-07,-2.14042e-10,1.53165e-14,4637.88,-91.6007], Tmin=(903.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C=C(C)C(=C)O(4205)',
    structure = SMILES('[CH2]C=C(C)C(=C)O'),
    E0 = (-75.5322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404469,0.0666345,-2.9704e-05,-1.76246e-08,1.43488e-11,-8943.57,24.108], Tmin=(100,'K'), Tmax=(937.055,'K')), NASAPolynomial(coeffs=[18.4317,0.0207826,-6.09127e-06,9.95705e-10,-6.86413e-14,-13687.5,-68.9767], Tmin=(937.055,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][CH][C](C)C([CH2])[O](4206)',
    structure = SMILES('[CH2][CH][C](C)C([CH2])[O]'),
    E0 = (627.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,274.339,991.087,2160.91],'cm^-1')),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02249,0.0646934,-4.81139e-05,1.89014e-08,-3.11709e-12,75594,31.8349], Tmin=(100,'K'), Tmax=(1381.09,'K')), NASAPolynomial(coeffs=[11.6828,0.0338185,-1.45808e-05,2.71467e-09,-1.87042e-13,72649.4,-23.0371], Tmin=(1381.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])[O](4207)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[O]'),
    E0 = (680.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,248.528,985.032,2220.78],'cm^-1')),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752264,0.073098,-7.18049e-05,4.06589e-08,-9.63206e-12,81917.4,34.4019], Tmin=(100,'K'), Tmax=(1005.11,'K')), NASAPolynomial(coeffs=[10.3935,0.0347274,-1.4539e-05,2.67408e-09,-1.83712e-13,79979.3,-12.1603], Tmin=(1005.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C](C)[O](4208)',
    structure = SMILES('[CH2][CH]C([CH2])[C](C)[O]'),
    E0 = (645.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,209.961,356.282,3404.73],'cm^-1')),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00308,0.0727402,-6.98797e-05,2.73547e-08,3.97349e-12,77698.5,32.3443], Tmin=(100,'K'), Tmax=(602.304,'K')), NASAPolynomial(coeffs=[7.94194,0.0390944,-1.7059e-05,3.17111e-09,-2.17961e-13,76637.1,0.513732], Tmin=(602.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[C]1CO1(4209)',
    structure = SMILES('[CH2][CH]C(C)[C]1CO1'),
    E0 = (399.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20291,0.0612063,-4.15686e-05,8.4838e-09,3.94376e-12,48141.3,28.7243], Tmin=(100,'K'), Tmax=(773.635,'K')), NASAPolynomial(coeffs=[8.62059,0.0343241,-1.16862e-05,1.89758e-09,-1.20875e-13,46650.4,-7.37634], Tmin=(773.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1C[C]([O])C1C(4210)',
    structure = SMILES('[CH2]C1C[C]([O])C1C'),
    E0 = (387.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20805,0.0473345,1.27246e-05,-5.0522e-08,2.30267e-11,46758.4,25.8006], Tmin=(100,'K'), Tmax=(962.211,'K')), NASAPolynomial(coeffs=[14.2405,0.0273195,-9.32975e-06,1.65666e-09,-1.17118e-13,42668.9,-44.7894], Tmin=(962.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl)"""),
)

species(
    label = 'CC1[CH]CC[C]1[O](4211)',
    structure = SMILES('CC1[CH]CC[C]1[O]'),
    E0 = (304.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68524,0.0360889,3.5044e-05,-6.59622e-08,2.62965e-11,36727.1,26.0286], Tmin=(100,'K'), Tmax=(990.611,'K')), NASAPolynomial(coeffs=[11.8597,0.0309918,-1.17288e-05,2.18684e-09,-1.57035e-13,32945.7,-31.8732], Tmin=(990.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1C(C)C1([CH2])[O](4212)',
    structure = SMILES('[CH2]C1C(C)C1([CH2])[O]'),
    E0 = (408.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596307,0.0626038,-2.25276e-05,-2.16574e-08,1.49012e-11,49229.9,25.7244], Tmin=(100,'K'), Tmax=(944.379,'K')), NASAPolynomial(coeffs=[17.1068,0.022775,-7.07934e-06,1.19017e-09,-8.23406e-14,44769.2,-60.0912], Tmin=(944.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(Isobutyl) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])C[CH]C1C(4213)',
    structure = SMILES('[CH2]C1([O])C[CH]C1C'),
    E0 = (395.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943387,0.0558524,-1.48782e-05,-1.78456e-08,9.93454e-12,47698,26.0339], Tmin=(100,'K'), Tmax=(1040.65,'K')), NASAPolynomial(coeffs=[13.8213,0.0301014,-1.19921e-05,2.2353e-09,-1.579e-13,43731.8,-42.7862], Tmin=(1040.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(cyclobutane) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C=C(C)C(=C)[O](4214)',
    structure = SMILES('[CH2]C=C(C)C(=C)[O]'),
    E0 = (62.2727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,365.277,365.38],'cm^-1')),
        HinderedRotor(inertia=(0.16107,'amu*angstrom^2'), symmetry=1, barrier=(15.0683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160691,'amu*angstrom^2'), symmetry=1, barrier=(15.0509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853478,'amu*angstrom^2'), symmetry=1, barrier=(80.7926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753129,0.0628798,-3.80576e-05,6.8208e-10,5.72992e-12,7614.26,24.357], Tmin=(100,'K'), Tmax=(956.437,'K')), NASAPolynomial(coeffs=[14.7682,0.0240139,-8.07407e-06,1.37016e-09,-9.26407e-14,4030.13,-47.3551], Tmin=(956.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.2727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([532.913,1391.12,1391.12,2779.21,3448.45,3448.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.948,0.0008276,8.34932e-06,-9.82634e-09,3.80104e-12,16425.4,0.336655], Tmin=(100,'K'), Tmax=(660.467,'K')), NASAPolynomial(coeffs=[3.2217,0.00522646,-1.64125e-06,2.58225e-10,-1.62579e-14,16521.3,3.53938], Tmin=(660.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C=CC(=C)[O](2637)',
    structure = SMILES('[CH2]C=CC(=C)[O]'),
    E0 = (99.8633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,298.677,302.296],'cm^-1')),
        HinderedRotor(inertia=(0.343709,'amu*angstrom^2'), symmetry=1, barrier=(22.5192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15335,'amu*angstrom^2'), symmetry=1, barrier=(74.3994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60392,0.0427776,-7.06824e-06,-2.71089e-08,1.55232e-11,12106.2,20.7036], Tmin=(100,'K'), Tmax=(925.173,'K')), NASAPolynomial(coeffs=[13.8391,0.015431,-4.15959e-06,6.4829e-10,-4.42362e-14,8748.69,-43.2829], Tmin=(925.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.8633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
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
    label = '[CH2][CH][CH]C(3856)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (450.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2810.32],'cm^-1')),
        HinderedRotor(inertia=(0.00169393,'amu*angstrom^2'), symmetry=1, barrier=(9.48541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010072,'amu*angstrom^2'), symmetry=1, barrier=(56.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100766,'amu*angstrom^2'), symmetry=1, barrier=(56.4312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79546,0.0335488,-4.51222e-05,5.0661e-08,-2.13657e-11,54172.4,20.3345], Tmin=(100,'K'), Tmax=(854.972,'K')), NASAPolynomial(coeffs=[-1.26875,0.0345785,-1.53757e-05,2.86235e-09,-1.94746e-13,55524.7,43.1492], Tmin=(854.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C(=C)[O](2638)',
    structure = SMILES('[CH2][CH]C=C([CH2])[O]'),
    E0 = (375.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,413.454,1819.17],'cm^-1')),
        HinderedRotor(inertia=(0.000978086,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000972178,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250014,'amu*angstrom^2'), symmetry=1, barrier=(30.5047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53573,0.0470622,-2.59979e-05,-3.43734e-09,5.97661e-12,45314.5,24.3495], Tmin=(100,'K'), Tmax=(952.509,'K')), NASAPolynomial(coeffs=[12.6541,0.0172656,-5.67952e-06,9.62703e-10,-6.55921e-14,42430,-32.7722], Tmin=(952.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][CH][C](C)C(=C)[O](4215)',
    structure = SMILES('[CH2][CH]C(C)=C([CH2])[O]'),
    E0 = (336.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.012328,'amu*angstrom^2'), symmetry=1, barrier=(3.41321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000393773,'amu*angstrom^2'), symmetry=1, barrier=(3.41345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0801695,'amu*angstrom^2'), symmetry=1, barrier=(22.1915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303417,'amu*angstrom^2'), symmetry=1, barrier=(84.0213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.575814,0.0683989,-6.10052e-05,2.90334e-08,-5.54211e-12,40651.1,28.5533], Tmin=(100,'K'), Tmax=(1264.99,'K')), NASAPolynomial(coeffs=[14.5169,0.0243159,-8.73228e-06,1.48477e-09,-9.76544e-14,37124.1,-41.9814], Tmin=(1264.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(=C)[O](4216)',
    structure = SMILES('[CH2][CH]C([CH2])C(=C)[O]'),
    E0 = (456.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,350.492,1313.27,1313.36],'cm^-1')),
        HinderedRotor(inertia=(0.00133196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00122576,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0519993,'amu*angstrom^2'), symmetry=1, barrier=(4.79708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475557,'amu*angstrom^2'), symmetry=1, barrier=(4.86684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796865,0.065193,-5.74532e-05,2.77245e-08,-5.42681e-12,54971.4,33.5014], Tmin=(100,'K'), Tmax=(1227.73,'K')), NASAPolynomial(coeffs=[12.8885,0.0257971,-9.31944e-06,1.58696e-09,-1.0435e-13,52002.5,-27.3141], Tmin=(1227.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C([O])C(C)[CH][CH2](4217)',
    structure = SMILES('[CH]=C([O])C(C)[CH][CH2]'),
    E0 = (498.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,3692.47],'cm^-1')),
        HinderedRotor(inertia=(0.996441,'amu*angstrom^2'), symmetry=1, barrier=(22.9101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14136,'amu*angstrom^2'), symmetry=1, barrier=(3.25015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269878,'amu*angstrom^2'), symmetry=1, barrier=(22.9111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00236766,'amu*angstrom^2'), symmetry=1, barrier=(22.9218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.821392,0.0656465,-5.66904e-05,2.60864e-08,-4.87928e-12,60022.7,32.0357], Tmin=(100,'K'), Tmax=(1273.74,'K')), NASAPolynomial(coeffs=[13.3298,0.0263659,-1.04325e-05,1.87546e-09,-1.27369e-13,56836.2,-31.3367], Tmin=(1273.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C[C](C)C(=C)[O](4218)',
    structure = SMILES('[CH2]CC(C)=C([CH2])[O]'),
    E0 = (195.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475609,0.0724947,-6.61358e-05,3.27392e-08,-6.54907e-12,23681.2,29.3249], Tmin=(100,'K'), Tmax=(1204.07,'K')), NASAPolynomial(coeffs=[13.9942,0.0275842,-1.01865e-05,1.76082e-09,-1.16967e-13,20425.8,-38.405], Tmin=(1204.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]CC([CH2])C(=C)[O](2504)',
    structure = SMILES('[CH2]CC([CH2])C(=C)[O]'),
    E0 = (261.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,405.999,406.002,406.003],'cm^-1')),
        HinderedRotor(inertia=(0.006215,'amu*angstrom^2'), symmetry=1, barrier=(0.726941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00621715,'amu*angstrom^2'), symmetry=1, barrier=(0.727194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00621712,'amu*angstrom^2'), symmetry=1, barrier=(0.727284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0582175,'amu*angstrom^2'), symmetry=1, barrier=(6.80968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667782,0.0656988,-4.36726e-05,6.65012e-09,3.58983e-12,31580.5,31.408], Tmin=(100,'K'), Tmax=(953.375,'K')), NASAPolynomial(coeffs=[14.2976,0.0260055,-8.74214e-06,1.46918e-09,-9.8094e-14,28186.7,-37.8661], Tmin=(953.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])[C](C)[CH]C(2742)',
    structure = SMILES('[CH2]C([O])=C(C)[CH]C'),
    E0 = (131.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,819.675],'cm^-1')),
        HinderedRotor(inertia=(0.0082949,'amu*angstrom^2'), symmetry=1, barrier=(4.05512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60485,'amu*angstrom^2'), symmetry=1, barrier=(13.9067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174309,'amu*angstrom^2'), symmetry=1, barrier=(4.00772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0293364,'amu*angstrom^2'), symmetry=1, barrier=(13.967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.604356,0.0669232,-4.81416e-05,1.39616e-08,-5.16752e-14,15965.6,25.9694], Tmin=(100,'K'), Tmax=(1028.35,'K')), NASAPolynomial(coeffs=[14.1014,0.0276215,-1.00656e-05,1.75773e-09,-1.18886e-13,12491.8,-42.9164], Tmin=(1028.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][CH][C](C)C(=C)O(4219)',
    structure = SMILES('[CH2][CH]C(C)=C([CH2])O'),
    E0 = (199.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.329807,0.0709925,-4.88714e-05,6.3165e-09,4.7098e-12,24088.9,27.9333], Tmin=(100,'K'), Tmax=(964.54,'K')), NASAPolynomial(coeffs=[17.2924,0.0225401,-7.56696e-06,1.29972e-09,-8.91374e-14,19798.3,-58.5682], Tmin=(964.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_S) + radical(C=C(O)CJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C(O)C(C)[CH][CH2](4220)',
    structure = SMILES('[CH]=C(O)C(C)[CH][CH2]'),
    E0 = (360.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,545.73],'cm^-1')),
        HinderedRotor(inertia=(0.00507696,'amu*angstrom^2'), symmetry=1, barrier=(1.07297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.054251,'amu*angstrom^2'), symmetry=1, barrier=(11.4654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542511,'amu*angstrom^2'), symmetry=1, barrier=(11.4654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542511,'amu*angstrom^2'), symmetry=1, barrier=(11.4654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.054251,'amu*angstrom^2'), symmetry=1, barrier=(11.4654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362098,0.0707336,-5.31852e-05,1.44484e-08,6.86782e-13,43469.8,32.1824], Tmin=(100,'K'), Tmax=(1021.67,'K')), NASAPolynomial(coeffs=[16.7077,0.0235846,-8.6954e-06,1.55678e-09,-1.07861e-13,39250.6,-51.3289], Tmin=(1021.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(=C)O(4221)',
    structure = SMILES('[CH2][CH]C([CH2])C(=C)O'),
    E0 = (318.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,886.59],'cm^-1')),
        HinderedRotor(inertia=(0.0755603,'amu*angstrom^2'), symmetry=1, barrier=(1.73728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0755616,'amu*angstrom^2'), symmetry=1, barrier=(1.73731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00311458,'amu*angstrom^2'), symmetry=1, barrier=(1.73726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700329,'amu*angstrom^2'), symmetry=1, barrier=(16.102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288679,'amu*angstrom^2'), symmetry=1, barrier=(16.1019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.496011,0.0684083,-4.73716e-05,7.48712e-09,3.84798e-12,38411.6,33.0798], Tmin=(100,'K'), Tmax=(960.186,'K')), NASAPolynomial(coeffs=[15.9848,0.0235017,-7.86494e-06,1.33537e-09,-9.04256e-14,34532.8,-45.7245], Tmin=(960.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C([O])C(C)C[CH2](4222)',
    structure = SMILES('[CH]=C([O])C(C)C[CH2]'),
    E0 = (303.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,363.824,366.879],'cm^-1')),
        HinderedRotor(inertia=(0.00128672,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0012325,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104192,'amu*angstrom^2'), symmetry=1, barrier=(9.65958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101371,'amu*angstrom^2'), symmetry=1, barrier=(9.67987,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.522635,0.0681599,-4.99765e-05,1.42671e-08,1.40578e-13,36639.2,30.5507], Tmin=(100,'K'), Tmax=(1028.25,'K')), NASAPolynomial(coeffs=[15.062,0.0260176,-9.53172e-06,1.68092e-09,-1.14727e-13,32887,-43.7046], Tmin=(1028.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(C)[CH]C(2744)',
    structure = SMILES('[CH]=C([O])C(C)[CH]C'),
    E0 = (292.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00838956,'amu*angstrom^2'), symmetry=1, barrier=(5.05059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660118,'amu*angstrom^2'), symmetry=1, barrier=(15.1774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660145,'amu*angstrom^2'), symmetry=1, barrier=(15.178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252088,'amu*angstrom^2'), symmetry=1, barrier=(15.1788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598771,0.0670083,-5.31394e-05,2.21779e-08,-3.74828e-12,35348.2,31.0538], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[14.7529,0.0267181,-1.01315e-05,1.77399e-09,-1.18255e-13,31370.3,-42.0467], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C(C)C[C]=O(2721)',
    structure = SMILES('[CH2][CH]C(C)C[C]=O'),
    E0 = (279.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1904.05,1904.05],'cm^-1')),
        HinderedRotor(inertia=(0.0839111,'amu*angstrom^2'), symmetry=1, barrier=(8.91661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00112577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839118,'amu*angstrom^2'), symmetry=1, barrier=(8.91662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839121,'amu*angstrom^2'), symmetry=1, barrier=(8.91661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839115,'amu*angstrom^2'), symmetry=1, barrier=(8.91662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23781,0.0656711,-6.39225e-05,4.1548e-08,-1.23157e-11,33698,30.9379], Tmin=(100,'K'), Tmax=(784.457,'K')), NASAPolynomial(coeffs=[5.8029,0.0423938,-1.94138e-05,3.72327e-09,-2.61537e-13,32981.8,10.0222], Tmin=(784.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH]C(C)[C]=O(4223)',
    structure = SMILES('[CH2][CH]C(C)[C]=O'),
    E0 = (307.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,4000],'cm^-1')),
        HinderedRotor(inertia=(0.426949,'amu*angstrom^2'), symmetry=1, barrier=(26.9945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.20416e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433321,'amu*angstrom^2'), symmetry=1, barrier=(26.9906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115374,'amu*angstrom^2'), symmetry=1, barrier=(7.2497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65405,0.0474606,-3.27582e-05,1.17544e-08,-1.74008e-12,37034.4,27.6069], Tmin=(100,'K'), Tmax=(1549.25,'K')), NASAPolynomial(coeffs=[11.038,0.0232323,-9.30007e-06,1.65999e-09,-1.11152e-13,34126.8,-21.7731], Tmin=(1549.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1CC(=O)C1C(4224)',
    structure = SMILES('[CH2]C1CC(=O)C1C'),
    E0 = (26.7903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67031,0.0371845,3.20747e-05,-6.39039e-08,2.61588e-11,3318.59,24.2862], Tmin=(100,'K'), Tmax=(971.807,'K')), NASAPolynomial(coeffs=[11.5698,0.0307612,-1.09898e-05,1.98269e-09,-1.40227e-13,-226.275,-31.5292], Tmin=(971.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.7903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C(C)C(C)=O(4225)',
    structure = SMILES('C=C[C](C)C(C)=O'),
    E0 = (-65.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30115,0.0493766,-7.25914e-06,-1.88312e-08,8.7484e-12,-7732.2,24.291], Tmin=(100,'K'), Tmax=(1086.24,'K')), NASAPolynomial(coeffs=[11.1589,0.0335945,-1.37996e-05,2.57279e-09,-1.80083e-13,-11084.3,-29.6545], Tmin=(1086.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]([CH2])O(4226)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])O'),
    E0 = (626.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,2769.9],'cm^-1')),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266716,0.0879048,-0.000122509,9.87791e-08,-3.19003e-11,75468.1,35.755], Tmin=(100,'K'), Tmax=(856.653,'K')), NASAPolynomial(coeffs=[9.36335,0.0361122,-1.55054e-05,2.81002e-09,-1.88113e-13,74251.4,-4.72817], Tmin=(856.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]1OC([CH2])C1C(4227)',
    structure = SMILES('[CH2][C]1OC([CH2])C1C'),
    E0 = (399.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.210456,0.0813494,-9.33375e-05,6.08314e-08,-1.55253e-11,48215.3,25.5875], Tmin=(100,'K'), Tmax=(1083.31,'K')), NASAPolynomial(coeffs=[12.2297,0.0284042,-8.16725e-06,1.11914e-09,-6.08886e-14,46113.8,-31.0406], Tmin=(1083.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1OC[CH]C1C(4228)',
    structure = SMILES('[CH2][C]1OC[CH]C1C'),
    E0 = (323.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34844,0.0507194,-5.32297e-06,-3.23735e-08,1.87862e-11,38987.9,23.0985], Tmin=(100,'K'), Tmax=(862.961,'K')), NASAPolynomial(coeffs=[11.1998,0.0285878,-7.75666e-06,1.10539e-09,-6.68011e-14,36411.4,-28.0536], Tmin=(862.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOCs) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH][C](C)C(C)=O(4229)',
    structure = SMILES('[CH2][CH]C(C)=C(C)[O]'),
    E0 = (178.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756858,0.0661289,-5.16275e-05,2.15722e-08,-3.70604e-12,21529.8,27.2146], Tmin=(100,'K'), Tmax=(1366.69,'K')), NASAPolynomial(coeffs=[13.1984,0.0297148,-1.16609e-05,2.07629e-09,-1.39733e-13,18129.1,-36.6948], Tmin=(1366.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)=O(4230)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)=O'),
    E0 = (304.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.437306,'amu*angstrom^2'), symmetry=1, barrier=(10.0545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182105,'amu*angstrom^2'), symmetry=1, barrier=(10.0549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00511079,'amu*angstrom^2'), symmetry=1, barrier=(33.0519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0499496,'amu*angstrom^2'), symmetry=1, barrier=(33.0523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43754,'amu*angstrom^2'), symmetry=1, barrier=(33.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02005,0.0696661,-6.71767e-05,3.87794e-08,-9.69296e-12,36683.9,30.6792], Tmin=(100,'K'), Tmax=(937.74,'K')), NASAPolynomial(coeffs=[8.33828,0.0384505,-1.72459e-05,3.28308e-09,-2.29969e-13,35311.4,-4.15681], Tmin=(937.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(CJC(C)C=O) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C=C([CH2])OC(4231)',
    structure = SMILES('[CH2][CH]C=C([CH2])OC'),
    E0 = (256.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611699,0.0632478,-2.97991e-05,-9.96846e-09,9.37452e-12,31024.4,27.8], Tmin=(100,'K'), Tmax=(984.023,'K')), NASAPolynomial(coeffs=[16.298,0.0244621,-8.75133e-06,1.56767e-09,-1.09996e-13,26728,-53.7693], Tmin=(984.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=C(O)CJ) + radical(RCCJ)"""),
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
    E0 = (250.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (250.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (701.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (445.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (393.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (665.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1084.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (716.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (705.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (259.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (314.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (649.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (743.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (670.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (482.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (387.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (367.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (408.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (395.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (281.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (270.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (319.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (435.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (512.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (611.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (548.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (667.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (709.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (387.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (379.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (359.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (408.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (552.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (401.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (347.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (325.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (524.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (744.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (259.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (314.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (651.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (493.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (323.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (368.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (386.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (570.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['ketene(1375)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['C=CC(C)C(=C)[O](2741)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]CC([CH2])=O(1905)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C(C)[CH]C(=C)[O](4200)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH2](502)', 'C=C([O])[CH]C(2381)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH2][CH]C(C)[C]=C(4201)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C]C(C)C(=C)[O](4202)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH][CH]C(C)C(=C)[O](4203)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C1OC(=C)C1C(4204)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C=C(C)C(=C)O(4205)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][C](C)C([CH2])[O](4206)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])[C](C)[O](4208)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2][CH]C(C)[C]1CO1(4209)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C1C[C]([O])C1C(4210)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(136.849,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 135.1 to 136.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['CC1[CH]CC[C]1[O](4211)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.49459e+10,'s^-1'), n=0.314867, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C1C(C)C1([CH2])[O](4212)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(157.645,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C1([O])C[CH]C1C(4213)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.36716e+08,'s^-1'), n=1.01412, Ea=(144.607,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 143.9 to 144.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]C=C(C)C(=C)[O](4214)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', '[CH2]C=CC(=C)[O](2637)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.139979,'m^3/(mol*s)'), n=2.09962, Ea=(33.817,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]=O(1376)', 'm1_allyl(186)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['ketene(1375)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH3(17)', '[CH2][CH][CH]C(=C)[O](2638)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=O(1376)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2][CH][C](C)C(=C)[O](4215)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C(=C)[O](4216)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH]=C([O])C(C)[CH][CH2](4217)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C[C](C)C(=C)[O](4218)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.23078e+06,'s^-1'), n=1.82328, Ea=(136.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]CC([CH2])C(=C)[O](2504)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['C=C([O])[C](C)[CH]C(2742)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.14627,'s^-1'), n=3.49598, Ea=(108.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_noH] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2][CH][C](C)C(=C)O(4219)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(O)C(C)[CH][CH2](4220)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]C([CH2])C(=C)O(4221)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.00963743,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C([O])C(C)C[CH2](4222)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([O])C(C)[CH]C(2744)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]C(C)C[C]=O(2721)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(T)(20)', '[CH2][CH]C(C)[C]=O(4223)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C1CC(=O)C1C(4224)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C=C(C)C(C)=O(4225)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])[C]([CH2])O(4226)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2][C]1OC([CH2])C1C(4227)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2][C]1OC[CH]C1C(4228)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.19156e+09,'s^-1'), n=0.640131, Ea=(72.3284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 71.1 to 72.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2][CH][C](C)C(C)=O(4229)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH]C([CH2])C(C)=O(4230)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]C=C([CH2])OC(4231)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = 'PDepNetwork #1205',
    isomers = [
        '[CH2][CH]C(C)C(=C)[O](2743)',
    ],
    reactants = [
        ('ketene(1375)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1205',
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


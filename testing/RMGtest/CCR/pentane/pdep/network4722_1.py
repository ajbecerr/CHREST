species(
    label = 'C=C([O])C(=C)[CH]O[O](19816)',
    structure = SMILES('C=C([O])C(=C)[CH]O[O]'),
    E0 = (167.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,365.271,365.782,366.4],'cm^-1')),
        HinderedRotor(inertia=(0.213009,'amu*angstrom^2'), symmetry=1, barrier=(20.2587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212055,'amu*angstrom^2'), symmetry=1, barrier=(20.2522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212784,'amu*angstrom^2'), symmetry=1, barrier=(20.2501,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0353323,0.0777402,-8.35661e-05,4.40777e-08,-8.89542e-12,20349.8,29.6976], Tmin=(100,'K'), Tmax=(1286.73,'K')), NASAPolynomial(coeffs=[19.9352,0.012289,-3.33858e-06,4.75854e-10,-2.85611e-14,15489.4,-70.5996], Tmin=(1286.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO)"""),
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
    label = 'C=[C]C(O[O])C(=C)[O](21044)',
    structure = SMILES('C=[C]C(O[O])C(=C)[O]'),
    E0 = (295.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,333.044,333.047,333.049],'cm^-1')),
        HinderedRotor(inertia=(0.12092,'amu*angstrom^2'), symmetry=1, barrier=(9.51774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120921,'amu*angstrom^2'), symmetry=1, barrier=(9.51774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120921,'amu*angstrom^2'), symmetry=1, barrier=(9.51774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4370.21,'J/mol'), sigma=(6.9753,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.62 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.713308,0.0746137,-8.98135e-05,5.72591e-08,-1.46184e-11,35658.6,31.2599], Tmin=(100,'K'), Tmax=(953.426,'K')), NASAPolynomial(coeffs=[12.6741,0.0244335,-1.08659e-05,2.05615e-09,-1.43488e-13,33377.8,-25.8733], Tmin=(953.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_S)"""),
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
    label = 'C=C([O])C(=C)[CH][O](22500)',
    structure = SMILES('[CH2]C(=C[O])C(=C)[O]'),
    E0 = (64.4117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30147,'amu*angstrom^2'), symmetry=1, barrier=(29.9235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29982,'amu*angstrom^2'), symmetry=1, barrier=(29.8854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.70548,0.0763064,-8.01811e-05,3.90501e-08,-6.92726e-12,7939.82,27.8012], Tmin=(100,'K'), Tmax=(1656,'K')), NASAPolynomial(coeffs=[21.9763,0.00356412,1.97243e-06,-5.70396e-10,4.25168e-14,2889.6,-85.6326], Tmin=(1656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH]C(=C)C(=C)[O](15787)',
    structure = SMILES('[CH]C(=C)C(=C)[O]'),
    E0 = (350.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.16229,'amu*angstrom^2'), symmetry=1, barrier=(49.7154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15557,'amu*angstrom^2'), symmetry=1, barrier=(49.5608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.45,'J/mol'), sigma=(6.08582,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.10 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1577,0.055669,-3.86782e-05,5.33014e-09,3.75373e-12,42315.3,20.6503], Tmin=(100,'K'), Tmax=(941.057,'K')), NASAPolynomial(coeffs=[13.7286,0.0187803,-6.2503e-06,1.03926e-09,-6.93381e-14,39216.8,-43.1255], Tmin=(941.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]O[O](21387)',
    structure = SMILES('[CH]O[O]'),
    E0 = (471.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.569664,'amu*angstrom^2'), symmetry=1, barrier=(13.0977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33696,0.0153854,-2.48163e-05,1.94592e-08,-5.80272e-12,56680.3,11.4536], Tmin=(100,'K'), Tmax=(834.588,'K')), NASAPolynomial(coeffs=[6.14825,0.00191155,-5.99818e-07,1.15174e-10,-8.22083e-15,56211.1,-1.6009], Tmin=(834.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CH2_triplet)"""),
)

species(
    label = 'C=[C]C(=C)[O](2700)',
    structure = SMILES('C=[C]C(=C)[O]'),
    E0 = (216.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44422,'amu*angstrom^2'), symmetry=1, barrier=(33.2055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6519,0.0455512,-4.93434e-05,2.74295e-08,-5.8038e-12,26168.2,16.7593], Tmin=(100,'K'), Tmax=(1302.68,'K')), NASAPolynomial(coeffs=[11.8537,0.00923612,-1.78234e-06,1.49166e-10,-4.0926e-15,23933.7,-33.5311], Tmin=(1302.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C(=C)[CH]O[O](22489)',
    structure = SMILES('C=[C]C(=C)[CH]O[O]'),
    E0 = (442.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.33696,'amu*angstrom^2'), symmetry=1, barrier=(53.7313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33174,'amu*angstrom^2'), symmetry=1, barrier=(30.6192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32318,'amu*angstrom^2'), symmetry=1, barrier=(30.4225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612806,0.0649171,-6.47387e-05,3.26741e-08,-6.3918e-12,53410.1,26.2253], Tmin=(100,'K'), Tmax=(1313.51,'K')), NASAPolynomial(coeffs=[16.4591,0.0143309,-4.30963e-06,6.53191e-10,-4.02402e-14,49448.2,-53.7803], Tmin=(1313.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJC=C)"""),
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
    label = 'C=C([O])C(=C)[C]O[O](22501)',
    structure = SMILES('[CH2]C(=[C]O[O])C(=C)[O]'),
    E0 = (439.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,276.061,276.113],'cm^-1')),
        HinderedRotor(inertia=(0.206334,'amu*angstrom^2'), symmetry=1, barrier=(11.1433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20634,'amu*angstrom^2'), symmetry=1, barrier=(11.142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206027,'amu*angstrom^2'), symmetry=1, barrier=(11.1417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390057,0.0813351,-0.000110986,7.45642e-08,-1.86513e-11,52951,30.3762], Tmin=(100,'K'), Tmax=(790.942,'K')), NASAPolynomial(coeffs=[14.9394,0.0169555,-6.33981e-06,1.06622e-09,-6.84296e-14,50361.7,-38.2226], Tmin=(790.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'C=C([O])C(=C)C1OO1(22502)',
    structure = SMILES('C=C([O])C(=C)C1OO1'),
    E0 = (14.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.704471,0.0781328,-8.03322e-05,3.94002e-08,-7.09608e-12,1922.53,29.6518], Tmin=(100,'K'), Tmax=(1622.63,'K')), NASAPolynomial(coeffs=[21.1214,0.00731578,8.60632e-07,-4.20184e-10,3.46805e-14,-2920.78,-79.3084], Tmin=(1622.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(dioxirane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OC(O[O])C1=C(21801)',
    structure = SMILES('C=C1OC(O[O])C1=C'),
    E0 = (59.5621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.638076,0.0635771,-4.49168e-05,6.7839e-09,3.3147e-12,7293.78,21.6134], Tmin=(100,'K'), Tmax=(1020.51,'K')), NASAPolynomial(coeffs=[17.575,0.0174088,-6.77293e-06,1.27864e-09,-9.23199e-14,2784.15,-65.5991], Tmin=(1020.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.5621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C([CH]O[O])C([CH2])[O](19807)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])[O]'),
    E0 = (648.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,895.683],'cm^-1')),
        HinderedRotor(inertia=(0.0956283,'amu*angstrom^2'), symmetry=1, barrier=(15.4523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672029,'amu*angstrom^2'), symmetry=1, barrier=(15.4513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210245,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271359,'amu*angstrom^2'), symmetry=1, barrier=(15.4522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202557,0.0774458,-8.31329e-05,4.38994e-08,-9.01416e-12,78118.2,34.4512], Tmin=(100,'K'), Tmax=(1194.96,'K')), NASAPolynomial(coeffs=[18.4095,0.0164997,-6.6286e-06,1.21751e-09,-8.4533e-14,73766.9,-56.6294], Tmin=(1194.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])[CH]O[O](10651)',
    structure = SMILES('[CH]=C([O])C([CH2])[CH]O[O]'),
    E0 = (579.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,240.979,2096.19],'cm^-1')),
        HinderedRotor(inertia=(0.218225,'amu*angstrom^2'), symmetry=1, barrier=(8.99285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00290301,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218232,'amu*angstrom^2'), symmetry=1, barrier=(8.9929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65144,'amu*angstrom^2'), symmetry=1, barrier=(68.0544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.135132,0.0887244,-0.000131163,1.00729e-07,-2.9926e-11,69823.2,34.0664], Tmin=(100,'K'), Tmax=(919.985,'K')), NASAPolynomial(coeffs=[13.3219,0.0224925,-8.66829e-06,1.4509e-09,-9.12874e-14,67773.4,-26.406], Tmin=(919.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])[C](C)[O](19808)',
    structure = SMILES('[CH]C([CH]O[O])=C(C)[O]'),
    E0 = (412.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416894,0.0755324,-7.41409e-05,3.85386e-08,-8.06197e-12,49761.3,31.1301], Tmin=(100,'K'), Tmax=(1153.81,'K')), NASAPolynomial(coeffs=[14.471,0.0268089,-1.07971e-05,1.93803e-09,-1.31441e-13,46518.2,-38.6833], Tmin=(1153.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])[C](C)[CH]O[O](22503)',
    structure = SMILES('[CH]C([O])=C(C)[CH]O[O]'),
    E0 = (412.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416894,0.0755324,-7.41409e-05,3.85386e-08,-8.06197e-12,49761.3,31.1301], Tmin=(100,'K'), Tmax=(1153.81,'K')), NASAPolynomial(coeffs=[14.471,0.0268089,-1.07971e-05,1.93803e-09,-1.31441e-13,46518.2,-38.6833], Tmin=(1153.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([O])[C]1CC1O[O](22504)',
    structure = SMILES('C=C([O])[C]1CC1O[O]'),
    E0 = (241.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07245,0.0509001,-9.1527e-06,-3.26454e-08,1.84448e-11,29125.2,27.7706], Tmin=(100,'K'), Tmax=(948.97,'K')), NASAPolynomial(coeffs=[17.5524,0.014611,-4.23126e-06,7.36748e-10,-5.4664e-14,24503.6,-58.7433], Tmin=(948.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJ(C)CO)"""),
)

species(
    label = 'C=C([CH]O[O])[C]1CO1(22505)',
    structure = SMILES('[CH2]C([CH]O[O])=C1CO1'),
    E0 = (288.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.62742,0.0572732,-1.10585e-05,-4.13621e-08,2.4015e-11,34820.2,26.2937], Tmin=(100,'K'), Tmax=(937.192,'K')), NASAPolynomial(coeffs=[22.0598,0.00844955,-1.17926e-06,1.69799e-10,-1.72669e-14,28929.8,-85.7079], Tmin=(937.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1[C]([O])CC1O[O](22506)',
    structure = SMILES('[CH2]C1=C([O])CC1O[O]'),
    E0 = (230.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790565,0.0569533,-2.06593e-05,-2.41125e-08,1.60804e-11,27894.1,27.51], Tmin=(100,'K'), Tmax=(951.783,'K')), NASAPolynomial(coeffs=[19.1481,0.0129982,-3.70124e-06,6.52696e-10,-4.94168e-14,22896.1,-68.0459], Tmin=(951.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1OC[C]1[CH]O[O](22507)',
    structure = SMILES('[CH2]C1=C([CH]O[O])CO1'),
    E0 = (270.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754013,0.0542306,-3.9879e-06,-4.71565e-08,2.56817e-11,32673.3,27.0195], Tmin=(100,'K'), Tmax=(938.486,'K')), NASAPolynomial(coeffs=[21.3836,0.00953262,-1.64018e-06,2.57677e-10,-2.35281e-14,26897.5,-81.3392], Tmin=(938.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(ROOJ) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])[C]1[CH]OOC1(22508)',
    structure = SMILES('[CH2]C([O])=C1[CH]OOC1'),
    E0 = (141.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47093,0.0339002,5.1034e-05,-9.95632e-08,4.29991e-11,17178.7,27.2694], Tmin=(100,'K'), Tmax=(941.317,'K')), NASAPolynomial(coeffs=[19.4174,0.0126031,-2.6147e-06,4.63087e-10,-4.09038e-14,11364.9,-71.1613], Tmin=(941.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C1[CH]OOC[C]1[O](22509)',
    structure = SMILES('[CH2]C1=C([O])COO[CH]1'),
    E0 = (121.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29464,0.0348337,5.66925e-05,-1.1144e-07,4.85421e-11,14726.4,22.7013], Tmin=(100,'K'), Tmax=(939.439,'K')), NASAPolynomial(coeffs=[22.2896,0.00804931,-5.08963e-07,9.42588e-11,-1.78063e-14,8018.97,-91.9802], Tmin=(939.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(36dihydro12dioxin) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([O])C(=C)C1O[O](21868)',
    structure = SMILES('[CH2]C1([O])C(=C)C1O[O]'),
    E0 = (458.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533223,0.0691926,-6.71328e-05,3.27342e-08,-6.27108e-12,55283.6,29.2153], Tmin=(100,'K'), Tmax=(1271.42,'K')), NASAPolynomial(coeffs=[16.6895,0.0183634,-7.16518e-06,1.29025e-09,-8.82217e-14,51175.3,-52.609], Tmin=(1271.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(ROOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1([CH]O[O])OC1=C(21827)',
    structure = SMILES('[CH2]C1([CH]O[O])OC1=C'),
    E0 = (395.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.532039,0.0893609,-0.000109048,6.3282e-08,-1.38053e-11,47790.7,27.4259], Tmin=(100,'K'), Tmax=(1229.96,'K')), NASAPolynomial(coeffs=[22.4743,0.00754039,-7.25328e-07,-5.90692e-11,9.96689e-15,42660.8,-86.1753], Tmin=(1229.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1(C(=C)[O])[CH]OO1(22510)',
    structure = SMILES('[CH2]C1(C(=C)[O])[CH]OO1'),
    E0 = (327.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.468748,0.0730061,-7.82318e-05,4.30856e-08,-9.31148e-12,39504.3,27.5855], Tmin=(100,'K'), Tmax=(1134.39,'K')), NASAPolynomial(coeffs=[15.7563,0.0191001,-6.95172e-06,1.19498e-09,-7.94532e-14,36036,-48.0955], Tmin=(1134.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=C(C)OJ) + radical(CCsJOO) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1([O])OO[CH]C1=C(22511)',
    structure = SMILES('[CH2]C1([O])OO[CH]C1=C'),
    E0 = (287.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975172,0.0393093,5.62393e-05,-1.20155e-07,5.4101e-11,34712.2,27.8527], Tmin=(100,'K'), Tmax=(928.263,'K')), NASAPolynomial(coeffs=[25.7627,0.00245073,2.76018e-06,-5.63645e-10,2.81767e-14,27096.5,-106.121], Tmin=(928.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CC(C)(O)OJ) + radical(C=CCJO) + radical(CJCOOH)"""),
)

species(
    label = 'C=C([O])C(=C)C=O(22512)',
    structure = SMILES('C=C([O])C(=C)C=O'),
    E0 = (-96.1204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2782.5,750,1395,475,1775,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.920077,'amu*angstrom^2'), symmetry=1, barrier=(21.1544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.923391,'amu*angstrom^2'), symmetry=1, barrier=(21.2306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510012,0.067351,-6.87242e-05,3.39519e-08,-6.46547e-12,-11427,22.5093], Tmin=(100,'K'), Tmax=(1293.48,'K')), NASAPolynomial(coeffs=[18.4591,0.0118447,-4.35554e-06,7.75924e-10,-5.32953e-14,-16070.3,-68.7034], Tmin=(1293.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.1204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = '[CH]=C([CH]O[O])C(=C)[O](19814)',
    structure = SMILES('[CH]=C([CH]O[O])C(=C)[O]'),
    E0 = (415.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,288.338,288.613],'cm^-1')),
        HinderedRotor(inertia=(0.348352,'amu*angstrom^2'), symmetry=1, barrier=(20.5584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348348,'amu*angstrom^2'), symmetry=1, barrier=(20.5578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348454,'amu*angstrom^2'), symmetry=1, barrier=(20.5554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211564,0.0775158,-8.73024e-05,4.5002e-08,-7.9313e-12,50055.6,29.2132], Tmin=(100,'K'), Tmax=(947.463,'K')), NASAPolynomial(coeffs=[18.9313,0.0115278,-3.48085e-06,5.51726e-10,-3.62707e-14,45922.9,-63.1776], Tmin=(947.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(=C)[CH]O[O](20860)',
    structure = SMILES('[CH]=C([O])C(=C)[CH]O[O]'),
    E0 = (415.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,288.338,288.613],'cm^-1')),
        HinderedRotor(inertia=(0.348352,'amu*angstrom^2'), symmetry=1, barrier=(20.5584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348348,'amu*angstrom^2'), symmetry=1, barrier=(20.5578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348454,'amu*angstrom^2'), symmetry=1, barrier=(20.5554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211564,0.0775158,-8.73024e-05,4.5002e-08,-7.9313e-12,50055.6,29.2132], Tmin=(100,'K'), Tmax=(947.463,'K')), NASAPolynomial(coeffs=[18.9313,0.0115278,-3.48085e-06,5.51726e-10,-3.62707e-14,45922.9,-63.1776], Tmin=(947.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(CO[O])C(=C)[O](19328)',
    structure = SMILES('[CH]=C(CO[O])C(=C)[O]'),
    E0 = (297.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.537864,'amu*angstrom^2'), symmetry=1, barrier=(12.3665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537531,'amu*angstrom^2'), symmetry=1, barrier=(12.3589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53785,'amu*angstrom^2'), symmetry=1, barrier=(12.3662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4369.1,'J/mol'), sigma=(6.96718,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.44 K, Pc=29.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.157089,0.0867217,-0.000120646,8.67813e-08,-2.44221e-11,35942.9,28.8147], Tmin=(100,'K'), Tmax=(876.674,'K')), NASAPolynomial(coeffs=[14.5225,0.0211739,-8.48787e-06,1.48718e-09,-9.78481e-14,33424.2,-38.5986], Tmin=(876.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C(=C)[CH]O[O](22513)',
    structure = SMILES('[CH]=C(O)C(=C)[CH]O[O]'),
    E0 = (277.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.0145,'amu*angstrom^2'), symmetry=1, barrier=(23.3253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01768,'amu*angstrom^2'), symmetry=1, barrier=(23.3984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01656,'amu*angstrom^2'), symmetry=1, barrier=(23.3728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01674,'amu*angstrom^2'), symmetry=1, barrier=(23.3767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.644467,0.0872441,-9.97952e-05,5.36737e-08,-1.08024e-11,33519.8,30.7849], Tmin=(100,'K'), Tmax=(1354.24,'K')), NASAPolynomial(coeffs=[23.8238,0.00618169,-2.7086e-07,-1.13587e-10,1.18985e-14,27698.8,-91.7036], Tmin=(1354.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C(=C)O(19818)',
    structure = SMILES('[CH]=C([CH]O[O])C(=C)O'),
    E0 = (277.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.0145,'amu*angstrom^2'), symmetry=1, barrier=(23.3253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01768,'amu*angstrom^2'), symmetry=1, barrier=(23.3984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01656,'amu*angstrom^2'), symmetry=1, barrier=(23.3728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01674,'amu*angstrom^2'), symmetry=1, barrier=(23.3767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.644467,0.0872441,-9.97952e-05,5.36737e-08,-1.08024e-11,33519.8,30.7849], Tmin=(100,'K'), Tmax=(1354.24,'K')), NASAPolynomial(coeffs=[23.8238,0.00618169,-2.7086e-07,-1.13587e-10,1.18985e-14,27698.8,-91.7036], Tmin=(1354.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(=C)CO[O](10648)',
    structure = SMILES('[CH]=C([O])C(=C)CO[O]'),
    E0 = (297.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.537864,'amu*angstrom^2'), symmetry=1, barrier=(12.3665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537531,'amu*angstrom^2'), symmetry=1, barrier=(12.3589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53785,'amu*angstrom^2'), symmetry=1, barrier=(12.3662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.157089,0.0867217,-0.000120646,8.67813e-08,-2.44221e-11,35942.9,28.8147], Tmin=(100,'K'), Tmax=(876.674,'K')), NASAPolynomial(coeffs=[14.5225,0.0211739,-8.48787e-06,1.48718e-09,-9.78481e-14,33424.2,-38.5986], Tmin=(876.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]OO)C(=C)[O](19815)',
    structure = SMILES('[CH]=C([CH]OO)C(=C)[O]'),
    E0 = (263.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,276.828,277.666],'cm^-1')),
        HinderedRotor(inertia=(0.499117,'amu*angstrom^2'), symmetry=1, barrier=(27.1379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493492,'amu*angstrom^2'), symmetry=1, barrier=(27.1595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502548,'amu*angstrom^2'), symmetry=1, barrier=(27.1454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494319,'amu*angstrom^2'), symmetry=1, barrier=(27.1488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.352791,0.084605,-9.45748e-05,5.07941e-08,-1.04028e-11,31798.2,30.4818], Tmin=(100,'K'), Tmax=(1244,'K')), NASAPolynomial(coeffs=[22.0156,0.0106911,-3.05089e-06,4.6006e-10,-2.90321e-14,26386.9,-81.6974], Tmin=(1244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(=C)[CH]OO(22514)',
    structure = SMILES('[CH]=C([O])C(=C)[CH]OO'),
    E0 = (263.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,276.828,277.666],'cm^-1')),
        HinderedRotor(inertia=(0.499117,'amu*angstrom^2'), symmetry=1, barrier=(27.1379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493492,'amu*angstrom^2'), symmetry=1, barrier=(27.1595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502548,'amu*angstrom^2'), symmetry=1, barrier=(27.1454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494319,'amu*angstrom^2'), symmetry=1, barrier=(27.1488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.352791,0.084605,-9.45748e-05,5.07941e-08,-1.04028e-11,31798.2,30.4818], Tmin=(100,'K'), Tmax=(1244,'K')), NASAPolynomial(coeffs=[22.0156,0.0106911,-3.05089e-06,4.6006e-10,-2.90321e-14,26386.9,-81.6974], Tmin=(1244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = 'C=C([CH]O[O])C[C]=O(19754)',
    structure = SMILES('C=C([CH]O[O])C[C]=O'),
    E0 = (215.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,180,1177.19],'cm^-1')),
        HinderedRotor(inertia=(0.0195733,'amu*angstrom^2'), symmetry=1, barrier=(19.0894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393197,'amu*angstrom^2'), symmetry=1, barrier=(19.112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.831423,'amu*angstrom^2'), symmetry=1, barrier=(19.1161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834897,'amu*angstrom^2'), symmetry=1, barrier=(19.1959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4149.48,'J/mol'), sigma=(6.68679,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=648.14 K, Pc=31.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731857,0.0627192,-5.10548e-05,1.96349e-08,-2.96228e-12,26036.1,30.8655], Tmin=(100,'K'), Tmax=(1582.77,'K')), NASAPolynomial(coeffs=[18.63,0.0174865,-8.18728e-06,1.57894e-09,-1.10307e-13,20370.4,-63.7011], Tmin=(1582.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(CCCJ=O)"""),
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
    label = 'C=C([C]=O)[CH]O[O](22515)',
    structure = SMILES('C=C([C]=O)[CH]O[O]'),
    E0 = (226.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.0157279,'amu*angstrom^2'), symmetry=1, barrier=(21.3781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926683,'amu*angstrom^2'), symmetry=1, barrier=(21.3063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926535,'amu*angstrom^2'), symmetry=1, barrier=(21.3029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83643,0.0517465,-5.21532e-05,2.72125e-08,-5.9075e-12,27327,24.0009], Tmin=(100,'K'), Tmax=(1080.28,'K')), NASAPolynomial(coeffs=[9.63885,0.0228567,-1.20394e-05,2.4578e-09,-1.78826e-13,25641.2,-14.2439], Tmin=(1080.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'C=C1C(=O)CC1O[O](21746)',
    structure = SMILES('C=C1C(=O)CC1O[O]'),
    E0 = (38.5411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42526,0.0489461,-2.52182e-05,2.48854e-09,1.02291e-12,4734.35,25.4449], Tmin=(100,'K'), Tmax=(1334.11,'K')), NASAPolynomial(coeffs=[13.2682,0.0251107,-1.15434e-05,2.21348e-09,-1.54534e-13,535.61,-38.9975], Tmin=(1334.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.5411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C([CH]O[O])[C]([CH2])O(19823)',
    structure = SMILES('[CH]C([CH]O[O])=C([CH2])O'),
    E0 = (433.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.483592,0.0859055,-9.02261e-05,4.69753e-08,-9.3914e-12,52340.9,33.551], Tmin=(100,'K'), Tmax=(1309.63,'K')), NASAPolynomial(coeffs=[21.0828,0.01549,-4.36873e-06,6.19493e-10,-3.64794e-14,47081.9,-74.8235], Tmin=(1309.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(ROOJ) + radical(C=CCJO) + radical(C=C(O)CJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]1OC(O[O])C1=C(22516)',
    structure = SMILES('[CH2]C1=C([CH2])C(O[O])O1'),
    E0 = (250.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.36026,0.0669906,-4.16997e-05,-6.78764e-09,1.0859e-11,30328.3,25.7756], Tmin=(100,'K'), Tmax=(954.364,'K')), NASAPolynomial(coeffs=[21.1072,0.0112498,-3.15199e-06,5.56919e-10,-4.26392e-14,24946.8,-80.7948], Tmin=(954.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(ROOJ) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[O]O[CH][C]1CCC1=O(22517)',
    structure = SMILES('[O]O[CH]C1=C([O])CC1'),
    E0 = (216.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05987,0.0493986,-6.76181e-07,-4.33693e-08,2.25662e-11,26138,28.0409], Tmin=(100,'K'), Tmax=(949.376,'K')), NASAPolynomial(coeffs=[18.4941,0.013661,-3.80527e-06,6.75774e-10,-5.20003e-14,21127.9,-64.1157], Tmin=(949.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C]1OOO[CH]C1=C(22518)',
    structure = SMILES('[CH2]C1=C([CH2])OOO[CH]1'),
    E0 = (384.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84924,0.054407,-1.0716e-05,-3.30946e-08,1.8523e-11,46320.8,24.954], Tmin=(100,'K'), Tmax=(967.447,'K')), NASAPolynomial(coeffs=[18.8761,0.0153457,-5.15137e-06,9.7048e-10,-7.35419e-14,41172.8,-69.9978], Tmin=(967.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(123trioxene) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C1([CH]O[O])CC1=O(21773)',
    structure = SMILES('[CH2]C1([CH]O[O])CC1=O'),
    E0 = (398.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.114142,0.0895737,-0.000130588,1.00561e-07,-3.04541e-11,48083.9,27.3064], Tmin=(100,'K'), Tmax=(859.776,'K')), NASAPolynomial(coeffs=[13.0985,0.0249928,-1.06372e-05,1.90676e-09,-1.26513e-13,46005.5,-32.4767], Tmin=(859.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)2C=O)"""),
)

species(
    label = '[CH]=C([CH]O[O])C(C)=O(19826)',
    structure = SMILES('[CH]=C([CH]O[O])C(C)=O'),
    E0 = (264.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,375,552.5,462.5,1710,3025,407.5,1350,352.5,3120,650,792.5,1650,2558.72],'cm^-1')),
        HinderedRotor(inertia=(0.149884,'amu*angstrom^2'), symmetry=1, barrier=(8.45671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488485,'amu*angstrom^2'), symmetry=1, barrier=(27.7423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493432,'amu*angstrom^2'), symmetry=1, barrier=(27.7062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00595114,'amu*angstrom^2'), symmetry=1, barrier=(27.6871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44217,0.0612582,-5.69957e-05,2.92819e-08,-6.46959e-12,31951.4,28.2339], Tmin=(100,'K'), Tmax=(1045.1,'K')), NASAPolynomial(coeffs=[8.88687,0.0327642,-1.60986e-05,3.19343e-09,-2.28844e-13,30395.4,-8.01075], Tmin=(1045.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C[C]=CO[O](21045)',
    structure = SMILES('C=C([O])C[C]=CO[O]'),
    E0 = (317.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,194.061,196.823,197.618],'cm^-1')),
        HinderedRotor(inertia=(0.276489,'amu*angstrom^2'), symmetry=1, barrier=(7.55233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278491,'amu*angstrom^2'), symmetry=1, barrier=(7.55641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272877,'amu*angstrom^2'), symmetry=1, barrier=(7.55291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4425.45,'J/mol'), sigma=(7.00962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.25 K, Pc=29.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528941,0.0746253,-8.84397e-05,5.39374e-08,-1.28791e-11,38318.3,31.9761], Tmin=(100,'K'), Tmax=(1029.93,'K')), NASAPolynomial(coeffs=[14.9298,0.0186956,-6.98313e-06,1.21107e-09,-8.06107e-14,35351.9,-37.9246], Tmin=(1029.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])[C]=CO[O](22519)',
    structure = SMILES('C=C([O])[C]=CO[O]'),
    E0 = (284.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.851508,'amu*angstrom^2'), symmetry=1, barrier=(19.5778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853301,'amu*angstrom^2'), symmetry=1, barrier=(19.6191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705573,0.0716595,-0.000102932,7.19251e-08,-1.89436e-11,34341,24.8172], Tmin=(100,'K'), Tmax=(1042.57,'K')), NASAPolynomial(coeffs=[14.8925,0.00937167,-2.01113e-06,1.62938e-10,-2.29854e-15,31809.9,-42.17], Tmin=(1042.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1OCC1=CO[O](21790)',
    structure = SMILES('C=C1OCC1=CO[O]'),
    E0 = (110.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813278,0.0579205,-2.69312e-05,-1.56965e-08,1.27376e-11,13437.8,23.4018], Tmin=(100,'K'), Tmax=(952.958,'K')), NASAPolynomial(coeffs=[18.2427,0.0141207,-4.20134e-06,7.3189e-10,-5.36583e-14,8782.75,-66.8398], Tmin=(952.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ)"""),
)

species(
    label = 'C=C([O])C1=COOC1(19179)',
    structure = SMILES('C=C([O])C1=COOC1'),
    E0 = (-20.2301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.800032,0.0562534,-1.7197e-05,-2.93505e-08,1.85405e-11,-2304.96,25.1631], Tmin=(100,'K'), Tmax=(939.788,'K')), NASAPolynomial(coeffs=[19.5314,0.011877,-2.78884e-06,4.52799e-10,-3.48536e-14,-7386.69,-72.3468], Tmin=(939.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.2301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1=COOOC1=C(22520)',
    structure = SMILES('[CH2]C1=COOOC1=C'),
    E0 = (218.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41726,0.0304789,7.16555e-05,-1.21848e-07,4.98145e-11,26380.3,20.4176], Tmin=(100,'K'), Tmax=(967.385,'K')), NASAPolynomial(coeffs=[21.2645,0.0141101,-4.83096e-06,1.06361e-09,-9.1333e-14,19466.3,-90.5643], Tmin=(967.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([O])C1=COC1(19805)',
    structure = SMILES('C=C([O])C1=COC1'),
    E0 = (-21.0617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.946056,0.0445569,2.65805e-05,-8.91491e-08,4.4213e-11,-2402.01,19.1826], Tmin=(100,'K'), Tmax=(909.222,'K')), NASAPolynomial(coeffs=[25.3213,-0.00269556,5.57837e-06,-1.19145e-09,7.74859e-14,-9313.85,-109.728], Tmin=(909.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.0617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1=COOC1=C(22521)',
    structure = SMILES('C=C1[CH]OOC1=C'),
    E0 = (177.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07508,0.0179646,8.41238e-05,-1.29226e-07,5.24059e-11,21411.8,19.7967], Tmin=(100,'K'), Tmax=(946.716,'K')), NASAPolynomial(coeffs=[18.4565,0.0105575,-2.06813e-06,4.29047e-10,-4.25445e-14,15540.4,-72.9656], Tmin=(946.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C([CH2])[O](22522)',
    structure = SMILES('[CH2]C(=[C]O[O])C([CH2])[O]'),
    E0 = (672.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,271.6,840.951],'cm^-1')),
        HinderedRotor(inertia=(0.0229541,'amu*angstrom^2'), symmetry=1, barrier=(11.516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500898,'amu*angstrom^2'), symmetry=1, barrier=(11.5166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00722642,'amu*angstrom^2'), symmetry=1, barrier=(11.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500867,'amu*angstrom^2'), symmetry=1, barrier=(11.5159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484422,0.0801008,-0.000103151,6.96189e-08,-1.866e-11,81009.3,35.2413], Tmin=(100,'K'), Tmax=(912.949,'K')), NASAPolynomial(coeffs=[13.3545,0.023714,-1.05097e-05,1.97209e-09,-1.36523e-13,78659.2,-25.678], Tmin=(912.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])[C](C)[O](22523)',
    structure = SMILES('[CH2]C([C]O[O])=C(C)[O]'),
    E0 = (537.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,218.463,219.749,2481.34],'cm^-1')),
        HinderedRotor(inertia=(0.419312,'amu*angstrom^2'), symmetry=1, barrier=(14.7403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434178,'amu*angstrom^2'), symmetry=1, barrier=(14.74,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43081,'amu*angstrom^2'), symmetry=1, barrier=(14.7412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.431659,'amu*angstrom^2'), symmetry=1, barrier=(14.7388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.14838,0.0882433,-0.000124264,9.04228e-08,-2.58748e-11,64739.8,29.4062], Tmin=(100,'K'), Tmax=(860.044,'K')), NASAPolynomial(coeffs=[14.233,0.0227331,-1.00016e-05,1.84647e-09,-1.25702e-13,62317.3,-36.4195], Tmin=(860.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C([O])[C]([CH2])CO[O](10650)',
    structure = SMILES('[CH]C([O])=C([CH2])CO[O]'),
    E0 = (446.847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39821,0.081761,-9.69995e-05,6.41682e-08,-1.72528e-11,53870.7,30.8878], Tmin=(100,'K'), Tmax=(903.199,'K')), NASAPolynomial(coeffs=[11.7669,0.0314101,-1.33741e-05,2.43972e-09,-1.65805e-13,51817.2,-22.8013], Tmin=(903.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C](C(=C)[O])C1OO1(22524)',
    structure = SMILES('[CH2]C([O])=C([CH2])C1OO1'),
    E0 = (198.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.646554,0.0783286,-8.116e-05,4.03275e-08,-7.37284e-12,24104.9,32.3561], Tmin=(100,'K'), Tmax=(1593.82,'K')), NASAPolynomial(coeffs=[20.8738,0.00783135,7.0469e-07,-4.05565e-10,3.4452e-14,19339.2,-74.9287], Tmin=(1593.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(dioxirane) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C1([O])CC1=CO[O](21847)',
    structure = SMILES('[CH2]C1([O])CC1=CO[O]'),
    E0 = (475.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493958,0.0698079,-7.06276e-05,3.62549e-08,-7.27478e-12,57330.3,30.6134], Tmin=(100,'K'), Tmax=(1221.96,'K')), NASAPolynomial(coeffs=[16.5623,0.0172098,-6.06213e-06,1.0301e-09,-6.82372e-14,53403.3,-50.1281], Tmin=(1221.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(ROOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(=[C]OO)C(=C)[O](22525)',
    structure = SMILES('[CH2]C(=[C]OO)C(=C)[O]'),
    E0 = (287.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,314.744,315.187],'cm^-1')),
        HinderedRotor(inertia=(0.220276,'amu*angstrom^2'), symmetry=1, barrier=(15.5402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220194,'amu*angstrom^2'), symmetry=1, barrier=(15.5387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220129,'amu*angstrom^2'), symmetry=1, barrier=(15.5423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55472,'amu*angstrom^2'), symmetry=1, barrier=(109.546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0751257,0.0856291,-0.000109279,7.00906e-08,-1.74748e-11,34682.7,30.7414], Tmin=(100,'K'), Tmax=(989.618,'K')), NASAPolynomial(coeffs=[16.7553,0.0182082,-7.08671e-06,1.24767e-09,-8.35374e-14,31381.3,-49.5569], Tmin=(989.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'C=C([O])C(C)=[C]O[O](22526)',
    structure = SMILES('C=C([O])C(C)=[C]O[O]'),
    E0 = (287.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301871,0.0850121,-0.00012444,9.56265e-08,-2.86767e-11,34731.2,30.4039], Tmin=(100,'K'), Tmax=(890.251,'K')), NASAPolynomial(coeffs=[12.8026,0.0227263,-9.18527e-06,1.59759e-09,-1.03709e-13,32747.9,-27.0902], Tmin=(890.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C(=C)O(22527)',
    structure = SMILES('[CH2]C(=[C]O[O])C(=C)O'),
    E0 = (301.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.706532,'amu*angstrom^2'), symmetry=1, barrier=(16.2446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.7064,'amu*angstrom^2'), symmetry=1, barrier=(16.2415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704961,'amu*angstrom^2'), symmetry=1, barrier=(16.2084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.705198,'amu*angstrom^2'), symmetry=1, barrier=(16.2139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.080296,0.0866634,-0.000108858,6.55386e-08,-1.46268e-11,36398.4,30.5552], Tmin=(100,'K'), Tmax=(911.265,'K')), NASAPolynomial(coeffs=[18.8625,0.0132622,-4.08184e-06,6.25227e-10,-3.88144e-14,32541.2,-61.2936], Tmin=(911.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])C(C)=CO[O](22528)',
    structure = SMILES('[CH]=C([O])C(C)=CO[O]'),
    E0 = (295.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.631647,'amu*angstrom^2'), symmetry=1, barrier=(14.5228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632953,'amu*angstrom^2'), symmetry=1, barrier=(14.5528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632791,'amu*angstrom^2'), symmetry=1, barrier=(14.5491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0304829,0.0865999,-0.000113416,7.28644e-08,-1.76157e-11,35629.5,28.6545], Tmin=(100,'K'), Tmax=(867.421,'K')), NASAPolynomial(coeffs=[17.2075,0.01573,-5.28425e-06,8.42641e-10,-5.27638e-14,32335.8,-53.5801], Tmin=(867.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[O]OC=C1CCC1=O(21734)',
    structure = SMILES('[O]OC=C1CCC1=O'),
    E0 = (55.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41331,0.0492758,-2.7872e-05,5.11828e-09,3.23868e-13,6779.86,26.7424], Tmin=(100,'K'), Tmax=(1321.79,'K')), NASAPolynomial(coeffs=[12.8319,0.0244284,-1.06908e-05,2.00897e-09,-1.38951e-13,2913.26,-34.7392], Tmin=(1321.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C1=COOCC1=O(22529)',
    structure = SMILES('[CH2]C1=COOCC1=O'),
    E0 = (-23.2753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54136,0.0356347,3.85735e-05,-7.12857e-08,2.73339e-11,-2694.94,19.3876], Tmin=(100,'K'), Tmax=(1045.11,'K')), NASAPolynomial(coeffs=[15.6379,0.026606,-1.29454e-05,2.70762e-09,-2.04785e-13,-8094.81,-60.9795], Tmin=(1045.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.2753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclohexane) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2]C1=COCC1=O(22530)',
    structure = SMILES('[CH2]C1=COCC1=O'),
    E0 = (-102.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92695,0.0172378,9.58934e-05,-1.50377e-07,6.23562e-11,-12238.7,18.2566], Tmin=(100,'K'), Tmax=(934.545,'K')), NASAPolynomial(coeffs=[22.1519,0.00319955,2.01443e-06,-3.64712e-10,1.1433e-14,-19186.2,-94.8928], Tmin=(934.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2][C](O)C([CH2])=[C]O[O](22531)',
    structure = SMILES('[CH2]C(O)=C([CH2])[C]O[O]'),
    E0 = (558.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,338.656,338.656],'cm^-1')),
        HinderedRotor(inertia=(0.171375,'amu*angstrom^2'), symmetry=1, barrier=(13.9474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30743,'amu*angstrom^2'), symmetry=1, barrier=(25.0202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171375,'amu*angstrom^2'), symmetry=1, barrier=(13.9474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171375,'amu*angstrom^2'), symmetry=1, barrier=(13.9474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171375,'amu*angstrom^2'), symmetry=1, barrier=(13.9474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425453,0.0948737,-0.000127724,8.30283e-08,-2.06318e-11,67305,30.6476], Tmin=(100,'K'), Tmax=(1000.05,'K')), NASAPolynomial(coeffs=[20.1675,0.0125048,-4.17563e-06,6.6546e-10,-4.18427e-14,63186.3,-68.7025], Tmin=(1000.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(ROOJ) + radical(Allyl_P) + radical(C=C(O)CJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C(C)=O(22532)',
    structure = SMILES('[CH2]C(=[C]O[O])C(C)=O'),
    E0 = (293.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,1685,370,265.242],'cm^-1')),
        HinderedRotor(inertia=(0.122638,'amu*angstrom^2'), symmetry=1, barrier=(6.09488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00159774,'amu*angstrom^2'), symmetry=1, barrier=(6.09692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122573,'amu*angstrom^2'), symmetry=1, barrier=(6.09815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121401,'amu*angstrom^2'), symmetry=1, barrier=(6.10403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998324,0.0729033,-0.000107645,9.49383e-08,-3.38908e-11,35457.3,30.826], Tmin=(100,'K'), Tmax=(801.774,'K')), NASAPolynomial(coeffs=[6.50614,0.0349446,-1.70222e-05,3.28282e-09,-2.28309e-13,34911,7.57155], Tmin=(801.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=C(C=O)CJ) + radical(C=CJO)"""),
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
    E0 = (189.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (389.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (583.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (397.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (743.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (962.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (650.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (175.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (176.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (711.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (642.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (437.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (437.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (363.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (354.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (293.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (294.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (304.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (231.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (458.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (397.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (327.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (287.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (167.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (435.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (422.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (598.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (626.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (626.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (442.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (469.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (321.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (342.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (296.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (414.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (461.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (664.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (176.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (458.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (293.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (298.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (384.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (400.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (309.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (487.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (722.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (176.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (175.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (218.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (250.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (420.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (735.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (545.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (471.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (214.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (475.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (429.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (423.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (345.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (339.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (176.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (176.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (265.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (566.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (338.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['ketene(1375)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(21.5376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 21.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(O[O])C(=C)[O](21044)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', 'C=C([O])C(=C)[CH][O](22500)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O2(2)', '[CH]C(=C)C(=C)[O](15787)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]O[O](21387)', 'C=[C]C(=C)[O](2700)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', 'C=[C]C(=C)[CH]O[O](22489)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C=C([O])C(=C)[C]O[O](22501)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C([O])C(=C)C1OO1(22502)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C1OC(O[O])C1=C(21801)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C([CH]O[O])C([CH2])[O](19807)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C([O])C([CH2])[CH]O[O](10651)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C([CH]O[O])[C](C)[O](19808)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])[C](C)[CH]O[O](22503)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C([O])[C]1CC1O[O](22504)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_csHO]
Euclidian distance = 2.8284271247461903
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C([CH]O[O])[C]1CO1(22505)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C1[C]([O])CC1O[O](22506)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.6836e+08,'s^-1'), n=0.948854, Ea=(125.32,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C1OC[C]1[CH]O[O](22507)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C([O])[C]1[CH]OOC1(22508)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.20877e+14,'s^-1'), n=-0.684234, Ea=(136.856,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secDe_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C1[CH]OOC[C]1[O](22509)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1([O])C(=C)C1O[O](21868)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(290.655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 289.7 to 290.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1([CH]O[O])OC1=C(21827)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(229.643,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1(C(=C)[O])[CH]OO1(22510)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(159.456,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra_2H;radadd_intra_O] for rate rule [R5;doublebond_intra_2H_secDe;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 159.1 to 159.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1([O])OO[CH]C1=C(22511)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.79781e+09,'s^-1'), n=0.611011, Ea=(119.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H;radadd_intra_O] + [R6;doublebond_intra_2H;radadd_intra] for rate rule [R6;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 115.5 to 119.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', 'C=C([O])C(=C)C=O(22512)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(21.0272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 20.9 to 21.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=O(1376)', 'C=C=CO[O](16806)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00267131,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['ketene(1375)', '[CH2][C]=CO[O](16807)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=O(1376)', '[CH2][C]=CO[O](16807)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH]=C([CH]O[O])C(=C)[O](19814)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH]=C([O])C(=C)[CH]O[O](20860)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(CO[O])C(=C)[O](19328)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C(O)C(=C)[CH]O[O](22513)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C([CH]O[O])C(=C)O(19818)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([O])C(=C)CO[O](10648)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C([CH]OO)C(=C)[O](19815)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C([O])C(=C)[CH]OO(22514)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([CH]O[O])C[C]=O(19754)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(T)(20)', 'C=C([C]=O)[CH]O[O](22515)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C1C(=O)CC1O[O](21746)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C([CH]O[O])[C]([CH2])O(19823)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2][C]1OC(O[O])C1=C(22516)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.14898e+08,'s^-1'), n=0.935645, Ea=(125.333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[O]O[CH][C]1CCC1=O(22517)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2][C]1OOO[CH]C1=C(22518)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.41956e+10,'s^-1'), n=0.267163, Ea=(216.167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 213.7 to 216.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1([CH]O[O])CC1=O(21773)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(232.727,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_(CO)_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C([CH]O[O])C(C)=O(19826)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C[C]=CO[O](21045)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CH2(T)(20)', 'C=C([O])[C]=CO[O](22519)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C1OCC1=CO[O](21790)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C([O])C1=COOC1(19179)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1=COOOC1=C(22520)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(50.4561,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSDSS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 45.2 to 50.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['O(4)', 'C=C([O])C1=COC1(19805)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['O(4)', '[CH2]C1=COOC1=C(22521)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(252.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_SSD;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation
Ea raised from 250.9 to 252.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(=[C]O[O])C([CH2])[O](22522)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(=[C]O[O])[C](C)[O](22523)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]=C([O])[C]([CH2])CO[O](10650)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2][C](C(=C)[O])C1OO1(22524)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1([O])CC1=CO[O](21847)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(307.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 306.7 to 307.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C(=[C]OO)C(=C)[O](22525)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['C=C([O])C(C)=[C]O[O](22526)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C(=[C]O[O])C(=C)O(22527)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH]=C([O])C(C)=CO[O](22528)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[O]OC=C1CCC1=O(21734)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C1=COOCC1=O(22529)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSDSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction63',
    reactants = ['C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['O(4)', '[CH2]C1=COCC1=O(22530)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(9.11625e+09,'s^-1'), n=0.55, Ea=(97.6966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;C_pri_rad_intra;OOJ] + [R4OO;C_pri_rad_intra;OO] for rate rule [R4OO;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2][C](O)C([CH2])=[C]O[O](22531)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C(=[C]O[O])C(C)=O(22532)'],
    products = ['C=C([O])C(=C)[CH]O[O](19816)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4722',
    isomers = [
        'C=C([O])C(=C)[CH]O[O](19816)',
    ],
    reactants = [
        ('ketene(1375)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4722',
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


species(
    label = '[CH]=CC(=C)[CH]O[O](20187)',
    structure = SMILES('[CH]=CC(=C)[CH]O[O]'),
    E0 = (491.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.06454,'amu*angstrom^2'), symmetry=1, barrier=(24.476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06233,'amu*angstrom^2'), symmetry=1, barrier=(24.425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06179,'amu*angstrom^2'), symmetry=1, barrier=(24.4127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906261,0.0585355,-4.07978e-05,2.64102e-09,5.49004e-12,59185.1,25.9758], Tmin=(100,'K'), Tmax=(965.349,'K')), NASAPolynomial(coeffs=[17.0331,0.0134073,-4.38561e-06,7.74773e-10,-5.5528e-14,55060.7,-56.4937], Tmin=(965.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(1342)',
    structure = SMILES('C#C'),
    E0 = (218.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,647.284,647.285,3592.16],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08745,0.00578572,8.56332e-06,-1.72824e-08,7.83594e-12,26273.8,4.4608], Tmin=(100,'K'), Tmax=(897.945,'K')), NASAPolynomial(coeffs=[5.89068,0.00209055,4.88462e-08,-5.66903e-11,4.15085e-15,25415.9,-10.7352], Tmin=(897.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.303,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=CC([C]=C)O[O](21265)',
    structure = SMILES('[CH]=CC([C]=C)O[O]'),
    E0 = (619.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.545172,'amu*angstrom^2'), symmetry=1, barrier=(12.5346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5451,'amu*angstrom^2'), symmetry=1, barrier=(12.5329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545133,'amu*angstrom^2'), symmetry=1, barrier=(12.5337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3800.67,'J/mol'), sigma=(6.26596,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.66 K, Pc=35.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34284,0.0624999,-7.2969e-05,4.78649e-08,-1.30081e-11,74584.5,27.8499], Tmin=(100,'K'), Tmax=(884.089,'K')), NASAPolynomial(coeffs=[9.18151,0.0270335,-1.27931e-05,2.4868e-09,-1.75933e-13,73198.5,-9.00128], Tmin=(884.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=CC(=C)[CH][O](22477)',
    structure = SMILES('[CH]C=C([CH2])C=O'),
    E0 = (366.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,328.755,328.949,328.965],'cm^-1')),
        HinderedRotor(inertia=(0.648717,'amu*angstrom^2'), symmetry=1, barrier=(49.8482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648996,'amu*angstrom^2'), symmetry=1, barrier=(49.8479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64873,'amu*angstrom^2'), symmetry=1, barrier=(49.8488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69043,0.0484727,-3.11521e-05,9.43606e-09,-1.15024e-12,44143.6,20.5011], Tmin=(100,'K'), Tmax=(1841.21,'K')), NASAPolynomial(coeffs=[13.4331,0.0229623,-1.03696e-05,1.91121e-09,-1.28527e-13,39819.4,-43.3188], Tmin=(1841.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C=O)CJ) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C(=C)C=[CH](15985)',
    structure = SMILES('[CH]C(=C)C=[CH]'),
    E0 = (674.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09856,'amu*angstrom^2'), symmetry=1, barrier=(48.25,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09724,'amu*angstrom^2'), symmetry=1, barrier=(48.2197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3076.13,'J/mol'), sigma=(5.35439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.48 K, Pc=45.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78075,0.0401882,-8.76076e-06,-1.97193e-08,1.12783e-11,81164.5,18.0732], Tmin=(100,'K'), Tmax=(955.832,'K')), NASAPolynomial(coeffs=[12.0562,0.0178508,-6.13458e-06,1.06669e-09,-7.39864e-14,78256.3,-35.9736], Tmin=(955.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]=C[C]=C(13645)',
    structure = SMILES('[CH]C=C=C'),
    E0 = (516.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1515,'amu*angstrom^2'), symmetry=1, barrier=(49.4672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.021819,8.22328e-06,-2.14765e-08,8.5561e-12,62145.7,13.8367], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82081,0.0192338,-7.45619e-06,1.36535e-09,-9.5319e-14,60610.2,-9.33483], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=CC(=C)[C]O[O](22478)',
    structure = SMILES('[CH]=CC([CH2])=[C]O[O]'),
    E0 = (762.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.656713,'amu*angstrom^2'), symmetry=1, barrier=(15.0991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.656923,'amu*angstrom^2'), symmetry=1, barrier=(15.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657127,'amu*angstrom^2'), symmetry=1, barrier=(15.1086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96848,0.0664909,-8.38917e-05,5.4214e-08,-1.3696e-11,91802.1,27.9525], Tmin=(100,'K'), Tmax=(974.924,'K')), NASAPolynomial(coeffs=[13.2758,0.0159951,-6.19933e-06,1.08658e-09,-7.2413e-14,89402.4,-31.1107], Tmin=(974.924,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]=CC(=C)[CH]O[O](22479)',
    structure = SMILES('[C]=CC(=C)[CH]O[O]'),
    E0 = (802.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,314.349,314.454],'cm^-1')),
        HinderedRotor(inertia=(0.299703,'amu*angstrom^2'), symmetry=1, barrier=(21.1186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300322,'amu*angstrom^2'), symmetry=1, barrier=(21.1164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384105,'amu*angstrom^2'), symmetry=1, barrier=(27.0899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.829967,0.0626816,-6.55993e-05,3.37894e-08,-6.74077e-12,96590.7,26.2268], Tmin=(100,'K'), Tmax=(1233.58,'K')), NASAPolynomial(coeffs=[16.3151,0.0124695,-4.54288e-06,7.92566e-10,-5.3569e-14,92770.2,-51.7308], Tmin=(1233.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(802.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C1OO1(22480)',
    structure = SMILES('[CH]=CC(=C)C1OO1'),
    E0 = (337.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08922,0.0490439,-3.553e-06,-4.5726e-08,2.58892e-11,40720.2,22.8596], Tmin=(100,'K'), Tmax=(908.491,'K')), NASAPolynomial(coeffs=[19.8926,0.00634949,7.38619e-07,-2.96134e-10,1.96844e-14,35649,-75.1578], Tmin=(908.491,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cds_P)"""),
)

species(
    label = 'C=C1C=CC1O[O](22458)',
    structure = SMILES('C=C1C=CC1O[O]'),
    E0 = (265.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50828,0.0400683,1.02305e-05,-4.83741e-08,2.29582e-11,32023.4,20.6682], Tmin=(100,'K'), Tmax=(959.628,'K')), NASAPolynomial(coeffs=[16.6424,0.013097,-4.0571e-06,7.65953e-10,-5.96671e-14,27456,-60.3848], Tmin=(959.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = '[CH]=[C]C([CH2])[CH]O[O](22481)',
    structure = SMILES('[CH]=[C]C([CH2])[CH]O[O]'),
    E0 = (894.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.282612,'amu*angstrom^2'), symmetry=1, barrier=(6.49781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00454856,'amu*angstrom^2'), symmetry=1, barrier=(6.50438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282732,'amu*angstrom^2'), symmetry=1, barrier=(6.50056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.68878,'amu*angstrom^2'), symmetry=1, barrier=(61.8203,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772698,0.0781451,-0.000124397,1.0626e-07,-3.48966e-11,107634,30.6935], Tmin=(100,'K'), Tmax=(894.561,'K')), NASAPolynomial(coeffs=[8.56897,0.0271157,-1.1719e-05,2.08319e-09,-1.35994e-13,106886,-2.43437], Tmin=(894.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(894.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]1CC1O[O](22482)',
    structure = SMILES('[CH]=C[C]1CC1O[O]'),
    E0 = (565.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63525,0.0396364,4.44025e-06,-3.74745e-08,1.79765e-11,68053.9,24.5957], Tmin=(100,'K'), Tmax=(970.409,'K')), NASAPolynomial(coeffs=[14.4452,0.0165344,-5.75889e-06,1.07134e-09,-7.90359e-14,64169.3,-44.0254], Tmin=(970.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Cds_P)"""),
)

species(
    label = '[O]O[CH][C]1C=CC1(22483)',
    structure = SMILES('[O]O[CH]C1=C[CH]C1'),
    E0 = (457.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90251,0.0292932,3.88619e-05,-7.45039e-08,3.11729e-11,55130.3,22.352], Tmin=(100,'K'), Tmax=(963.857,'K')), NASAPolynomial(coeffs=[15.1828,0.0159437,-5.3577e-06,1.03593e-09,-8.03175e-14,50630.2,-51.2925], Tmin=(963.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(cyclobutene-allyl) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C[C]1[CH]OOC1(22484)',
    structure = SMILES('[CH]C=C1[CH]OOC1'),
    E0 = (430.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17014,0.0181703,8.60734e-05,-1.23184e-07,4.76101e-11,51854.7,23.6213], Tmin=(100,'K'), Tmax=(965.733,'K')), NASAPolynomial(coeffs=[14.3085,0.0234572,-8.44044e-06,1.6369e-09,-1.24999e-13,46919.2,-47.9307], Tmin=(965.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentane) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1([CH]O[O])C=C1(22485)',
    structure = SMILES('[CH2]C1([CH]O[O])C=C1'),
    E0 = (665.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2587,0.0631715,-6.82263e-05,3.98709e-08,-9.57306e-12,80083.6,23.4832], Tmin=(100,'K'), Tmax=(997.19,'K')), NASAPolynomial(coeffs=[10.5603,0.0258606,-1.21029e-05,2.35034e-09,-1.66594e-13,78228.5,-21.3657], Tmin=(997.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(ROOJ) + radical(CCsJOOH) + radical(Neopentyl)"""),
)

species(
    label = '[CH]=CC1([CH2])[CH]OO1(22486)',
    structure = SMILES('[CH]=CC1([CH2])[CH]OO1'),
    E0 = (651.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17563,0.0600276,-5.85599e-05,3.02361e-08,-6.28854e-12,78426.8,23.8947], Tmin=(100,'K'), Tmax=(1158.39,'K')), NASAPolynomial(coeffs=[12.3001,0.0216141,-8.81843e-06,1.60938e-09,-1.10425e-13,75849.5,-31.4102], Tmin=(1158.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = 'C#CC(=C)[CH]O[O](22487)',
    structure = SMILES('C#CC(=C)[CH]O[O]'),
    E0 = (417.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,410.035],'cm^-1')),
        HinderedRotor(inertia=(0.212133,'amu*angstrom^2'), symmetry=1, barrier=(25.2216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660149,'amu*angstrom^2'), symmetry=1, barrier=(78.5903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211709,'amu*angstrom^2'), symmetry=1, barrier=(25.2088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25943,0.0537122,-4.59926e-05,1.73536e-08,-1.75943e-12,50368.6,24.1727], Tmin=(100,'K'), Tmax=(1047.15,'K')), NASAPolynomial(coeffs=[13.8828,0.0156118,-5.91095e-06,1.0642e-09,-7.36899e-14,47170,-39.9584], Tmin=(1047.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=CC(=C)C=O(13714)',
    structure = SMILES('[CH]=CC(=C)C=O'),
    E0 = (227.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.06459,'amu*angstrom^2'), symmetry=1, barrier=(24.4771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06117,'amu*angstrom^2'), symmetry=1, barrier=(24.3984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30551,0.0499053,-3.22785e-05,9.60926e-10,4.21178e-12,27414.7,19.3092], Tmin=(100,'K'), Tmax=(1029.15,'K')), NASAPolynomial(coeffs=[15.6535,0.0127556,-5.26642e-06,1.03998e-09,-7.72081e-14,23475.5,-55.1139], Tmin=(1029.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(T)(1343)',
    structure = SMILES('[CH]=[CH]'),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([597.285,1197.51,1197.81,1199.36,2685.83,3758.62],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88059,-0.000843477,2.04095e-05,-2.51946e-08,9.47585e-12,71059.3,4.72114], Tmin=(100,'K'), Tmax=(935.365,'K')), NASAPolynomial(coeffs=[4.9214,0.00358273,-9.24444e-07,1.57274e-10,-1.19541e-14,70476.2,-2.3065], Tmin=(935.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=[C]C(=C)[CH]O[O](22488)',
    structure = SMILES('[CH]=C=C([CH2])[CH]O[O]'),
    E0 = (601.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.19253,'amu*angstrom^2'), symmetry=1, barrier=(64.5245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633345,'amu*angstrom^2'), symmetry=1, barrier=(34.2486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19048,'amu*angstrom^2'), symmetry=1, barrier=(64.4855,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793688,0.0626483,-6.55771e-05,3.45191e-08,-7.03206e-12,72440.9,27.0706], Tmin=(100,'K'), Tmax=(1242.29,'K')), NASAPolynomial(coeffs=[15.8385,0.0129064,-3.94719e-06,6.03723e-10,-3.74201e-14,68803.2,-48.3724], Tmin=(1242.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=CC(=[CH])[CH]O[O](20184)',
    structure = SMILES('[CH]=CC(=[CH])[CH]O[O]'),
    E0 = (738.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(1.11186,'amu*angstrom^2'), symmetry=1, barrier=(25.5639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11145,'amu*angstrom^2'), symmetry=1, barrier=(25.5544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11094,'amu*angstrom^2'), symmetry=1, barrier=(25.5426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50899,0.0658022,-7.02051e-05,3.60486e-08,-7.03643e-12,88919,27.8089], Tmin=(100,'K'), Tmax=(1349.9,'K')), NASAPolynomial(coeffs=[18.5264,0.0084969,-2.17592e-06,3.02222e-10,-1.81999e-14,84411.5,-63.1987], Tmin=(1349.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(=C)[CH]O[O](22489)',
    structure = SMILES('C=[C]C(=C)[CH]O[O]'),
    E0 = (442.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612806,0.0649171,-6.47387e-05,3.26741e-08,-6.3918e-12,53410.1,26.2253], Tmin=(100,'K'), Tmax=(1313.51,'K')), NASAPolynomial(coeffs=[16.4591,0.0143309,-4.30963e-06,6.53191e-10,-4.02402e-14,49448.2,-53.7803], Tmin=(1313.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]C(=C)CO[O](20188)',
    structure = SMILES('[CH]=C=C([CH2])CO[O]'),
    E0 = (483.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08458,0.0675651,-8.27274e-05,5.3124e-08,-1.23734e-11,58313.3,25.4461], Tmin=(100,'K'), Tmax=(718.761,'K')), NASAPolynomial(coeffs=[10.3347,0.0244063,-1.00205e-05,1.79062e-09,-1.1983e-13,56768.7,-17.6208], Tmin=(718.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=CC(=[CH])CO[O](19485)',
    structure = SMILES('[CH]=CC(=[CH])CO[O]'),
    E0 = (620.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(0.711795,'amu*angstrom^2'), symmetry=1, barrier=(16.3656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713106,'amu*angstrom^2'), symmetry=1, barrier=(16.3957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712245,'amu*angstrom^2'), symmetry=1, barrier=(16.3759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3799.01,'J/mol'), sigma=(6.25763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.40 K, Pc=35.18 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862228,0.070237,-8.70493e-05,5.6714e-08,-1.46435e-11,74788.6,25.9455], Tmin=(100,'K'), Tmax=(948.433,'K')), NASAPolynomial(coeffs=[12.6224,0.0206366,-8.60029e-06,1.56883e-09,-1.07014e-13,72557.9,-30.1675], Tmin=(948.433,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C=C(20189)',
    structure = SMILES('[CH]=C([CH]O[O])C=C'),
    E0 = (491.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.06219,'amu*angstrom^2'), symmetry=1, barrier=(24.4219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06246,'amu*angstrom^2'), symmetry=1, barrier=(24.428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06402,'amu*angstrom^2'), symmetry=1, barrier=(24.464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906261,0.0585355,-4.07978e-05,2.64102e-09,5.49004e-12,59185.1,25.9758], Tmin=(100,'K'), Tmax=(965.349,'K')), NASAPolynomial(coeffs=[17.0331,0.0134073,-4.38561e-06,7.74773e-10,-5.5528e-14,55060.7,-56.4937], Tmin=(965.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(=C)[CH]OO(22490)',
    structure = SMILES('[CH]=C=C([CH2])[CH]OO'),
    E0 = (449.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705992,0.0641747,-5.37161e-05,1.6039e-08,5.89791e-13,54162.8,26.626], Tmin=(100,'K'), Tmax=(982.891,'K')), NASAPolynomial(coeffs=[16.86,0.015493,-5.45667e-06,9.64184e-10,-6.7308e-14,50163.2,-55.2208], Tmin=(982.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=CC(=[CH])[CH]OO(20186)',
    structure = SMILES('[CH]=CC(=[CH])[CH]OO'),
    E0 = (586.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(1.31601,'amu*angstrom^2'), symmetry=1, barrier=(30.2576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31607,'amu*angstrom^2'), symmetry=1, barrier=(30.259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3142,'amu*angstrom^2'), symmetry=1, barrier=(30.2161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31578,'amu*angstrom^2'), symmetry=1, barrier=(30.2524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540866,0.0659678,-5.38077e-05,1.19898e-08,2.83869e-12,70635.6,26.9319], Tmin=(100,'K'), Tmax=(980.476,'K')), NASAPolynomial(coeffs=[19.1611,0.0117221,-4.04528e-06,7.46192e-10,-5.49182e-14,65940.4,-67.8565], Tmin=(980.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC[C]=CO[O](21266)',
    structure = SMILES('[CH]=CC[C]=CO[O]'),
    E0 = (641.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.499516,'amu*angstrom^2'), symmetry=1, barrier=(11.4848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499549,'amu*angstrom^2'), symmetry=1, barrier=(11.4856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499922,'amu*angstrom^2'), symmetry=1, barrier=(11.4942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3857.61,'J/mol'), sigma=(6.30207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.55 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21648,0.0618441,-6.93336e-05,4.16838e-08,-1.00659e-11,77241.7,28.3569], Tmin=(100,'K'), Tmax=(1005.71,'K')), NASAPolynomial(coeffs=[11.3987,0.021347,-8.93346e-06,1.64601e-09,-1.13336e-13,75193.6,-20.8243], Tmin=(1005.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=C[C]=CO[O](22491)',
    structure = SMILES('[CH]C=C=CO[O]'),
    E0 = (585.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.12634,'amu*angstrom^2'), symmetry=1, barrier=(48.8888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12731,'amu*angstrom^2'), symmetry=1, barrier=(48.911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72806,0.0498253,-4.99198e-05,2.77264e-08,-6.32435e-12,70452.5,21.9523], Tmin=(100,'K'), Tmax=(1051.92,'K')), NASAPolynomial(coeffs=[9.34205,0.0208728,-8.63484e-06,1.5618e-09,-1.0609e-13,68850.6,-15.1663], Tmin=(1051.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]OC=C1C=CC1(22444)',
    structure = SMILES('[O]OC=C1C=CC1'),
    E0 = (287.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54271,0.0373819,2.15759e-05,-6.53641e-08,3.08439e-11,34673.8,21.9938], Tmin=(100,'K'), Tmax=(930.098,'K')), NASAPolynomial(coeffs=[18.1414,0.00867055,-9.4088e-07,1.0374e-10,-1.20706e-14,29740.4,-66.8049], Tmin=(930.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = '[CH]=CC1=COOC1(19217)',
    structure = SMILES('[CH]=CC1=COOC1'),
    E0 = (302.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42736,0.0407212,1.29038e-05,-5.46405e-08,2.61679e-11,36544.1,22.5707], Tmin=(100,'K'), Tmax=(946.962,'K')), NASAPolynomial(coeffs=[17.8404,0.0109794,-2.69158e-06,4.84604e-10,-3.98655e-14,31660.6,-65.0901], Tmin=(946.962,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolene) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=COOC=C1(22492)',
    structure = SMILES('C=C1[CH]OOC=C1'),
    E0 = (155.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06154,0.0121899,0.000116021,-1.65098e-07,6.45607e-11,18854.3,16.2313], Tmin=(100,'K'), Tmax=(963.252,'K')), NASAPolynomial(coeffs=[20.6457,0.012088,-3.83734e-06,9.19879e-10,-8.51649e-14,11698.5,-91.2906], Tmin=(963.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=CC1=COC1(20181)',
    structure = SMILES('[CH]=CC1=COC1'),
    E0 = (302.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57006,0.0290743,5.64509e-05,-1.1404e-07,5.16164e-11,36447.2,16.6015], Tmin=(100,'K'), Tmax=(915.651,'K')), NASAPolynomial(coeffs=[23.6064,-0.00355411,5.65379e-06,-1.1546e-09,7.20625e-14,29743.9,-102.336], Tmin=(915.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=COC=C1(22493)',
    structure = SMILES('[CH2]C1=COC=C1'),
    E0 = (66.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.41943,0.0154755,6.72595e-05,-1.02257e-07,4.10085e-11,8097.82,16.5467], Tmin=(100,'K'), Tmax=(953.134,'K')), NASAPolynomial(coeffs=[15.0016,0.0114284,-3.10188e-06,6.26066e-10,-5.39111e-14,3484.66,-55.1689], Tmin=(953.134,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Furan) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]1C=CC1O[O](22494)',
    structure = SMILES('[CH2]C1=C[CH]C1O[O]'),
    E0 = (424.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21265,0.0460726,2.03825e-06,-4.27775e-08,2.12837e-11,51173.4,21.6363], Tmin=(100,'K'), Tmax=(968.73,'K')), NASAPolynomial(coeffs=[17.9148,0.0137619,-4.68755e-06,9.10025e-10,-7.05995e-14,46217.5,-67.2883], Tmin=(968.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C[C]([CH2])C1OO1(22495)',
    structure = SMILES('[CH]C=C([CH2])C1OO1'),
    E0 = (487.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13952,0.0499406,-2.45602e-06,-3.93065e-08,2.09576e-11,58733.1,24.7968], Tmin=(100,'K'), Tmax=(930.517,'K')), NASAPolynomial(coeffs=[16.0353,0.018851,-5.44332e-06,8.79922e-10,-6.09623e-14,54534.8,-53.6576], Tmin=(930.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(dioxirane) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC([CH2])=[C]OO(22496)',
    structure = SMILES('[CH]=CC([CH2])=[C]OO'),
    E0 = (610.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.818508,'amu*angstrom^2'), symmetry=1, barrier=(18.8191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54505,'amu*angstrom^2'), symmetry=1, barrier=(35.5237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.817419,'amu*angstrom^2'), symmetry=1, barrier=(18.7941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819335,'amu*angstrom^2'), symmetry=1, barrier=(18.8381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686871,0.0702973,-7.99353e-05,4.58138e-08,-1.02735e-11,73532.4,28.2036], Tmin=(100,'K'), Tmax=(1093.93,'K')), NASAPolynomial(coeffs=[15.2784,0.0169429,-6.77563e-06,1.22863e-09,-8.43112e-14,70340,-43.5021], Tmin=(1093.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(C)=[C]O[O](22497)',
    structure = SMILES('[CH]=CC(C)=[C]O[O]'),
    E0 = (610.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00133,0.0685783,-9.09501e-05,6.56106e-08,-1.88903e-11,73577.2,27.5561], Tmin=(100,'K'), Tmax=(851.03,'K')), NASAPolynomial(coeffs=[10.8504,0.0222861,-9.35734e-06,1.694e-09,-1.14141e-13,71900.8,-18.3712], Tmin=(851.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(C)=CO[O](22498)',
    structure = SMILES('[CH]=C=C(C)[CH]O[O]'),
    E0 = (449.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0428,0.0622237,-6.40168e-05,3.49306e-08,-7.59885e-12,54206.5,25.8958], Tmin=(100,'K'), Tmax=(1118.24,'K')), NASAPolynomial(coeffs=[12.7838,0.0202252,-7.67993e-06,1.34375e-09,-8.9947e-14,51580.7,-32.06], Tmin=(1118.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CCJO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C=C(22499)',
    structure = SMILES('[CH2]C(=[C]O[O])C=C'),
    E0 = (515.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,226.387],'cm^-1')),
        HinderedRotor(inertia=(0.38532,'amu*angstrom^2'), symmetry=1, barrier=(14.0132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385312,'amu*angstrom^2'), symmetry=1, barrier=(14.0132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385313,'amu*angstrom^2'), symmetry=1, barrier=(14.0132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0067,0.0634113,-6.88767e-05,3.90541e-08,-8.75276e-12,62083.9,27.4104], Tmin=(100,'K'), Tmax=(1090.71,'K')), NASAPolynomial(coeffs=[13.2762,0.0184138,-6.99234e-06,1.22801e-09,-8.24994e-14,59407.4,-32.848], Tmin=(1090.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
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
    E0 = (491.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (713.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (885.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (720.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1042.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (974.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1013.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (498.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (499.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (916.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (686.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (629.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (627.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (665.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (651.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (648.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (491.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (865.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (683.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1028.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (813.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (949.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (659.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (667.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (766.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (894.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (584.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (619.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (811.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1022.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (499.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (498.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (499.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (573.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (545.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (629.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (537.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (752.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (746.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (693.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (918.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['C2H2(1342)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC([C]=C)O[O](21265)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]=CC(=C)[CH][O](22477)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O2(2)', '[CH]C(=C)C=[CH](15985)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]O[O](21387)', '[CH]=C[C]=C(13645)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]=CC(=C)[C]O[O](22478)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[C]=CC(=C)[CH]O[O](22479)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=CC(=C)C1OO1(22480)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['C=C1C=CC1O[O](22458)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]C([CH2])[CH]O[O](22481)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=C[C]1CC1O[O](22482)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_csHO]
Euclidian distance = 2.8284271247461903
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[O]O[CH][C]1C=CC1(22483)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=C[C]1[CH]OOC1(22484)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.20877e+14,'s^-1'), n=-0.684234, Ea=(136.856,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secDe_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH2]C1([CH]O[O])C=C1(22485)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(174.884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=CC1([CH2])[CH]OO1(22486)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(160.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra_2H;radadd_intra_O] for rate rule [R5;doublebond_intra_2H_secDe;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C#CC(=C)[CH]O[O](22487)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.642e+08,'cm^3/(mol*s)'), n=1.548, Ea=(19.0205,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 195 used for Ct-Cd_Ct-H;HJ
Exact match found for rate rule [Ct-Cd_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O(4)', '[CH]=CC(=C)C=O(13714)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(21.0272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 20.9 to 21.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C2H2(T)(1343)', 'C=C=CO[O](16806)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00534261,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C2H2(1342)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H2(T)(1343)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.26444,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH]=[C]C(=C)[CH]O[O](22488)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH]=CC(=[CH])[CH]O[O](20184)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['C=[C]C(=C)[CH]O[O](22489)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=[C]C(=C)CO[O](20188)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00351592,'s^-1'), n=4.63833, Ea=(175.937,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeO;XH_out] + [R3H_SS_2Cd;C_rad_out_1H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=CC(=[CH])CO[O](19485)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C([CH]O[O])C=C(20189)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=[C]C(=C)[CH]OO(22490)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(36450.5,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_2;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=CC(=[CH])[CH]OO(20186)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC[C]=CO[O](21266)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(T)(20)', '[CH]=C[C]=CO[O](22491)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[O]OC=C1C=CC1(22444)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=CC1=COOC1(19217)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH2]C1=COOC=C1(22492)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['O(4)', '[CH]=CC1=COC1(20181)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['O(4)', '[CH2]C1=COC=C1(22493)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_DSD;Cd_pri_rad_in;OOJ]
Euclidian distance = 2.449489742783178
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH2][C]1C=CC1O[O](22494)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=C[C]([CH2])C1OO1(22495)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC([CH2])=[C]OO(22496)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=CC(C)=[C]O[O](22497)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=CC(=C)[CH]O[O](20187)'],
    products = ['[CH]=[C]C(C)=CO[O](22498)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(=[C]O[O])C=C(22499)'],
    products = ['[CH]=CC(=C)[CH]O[O](20187)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4979',
    isomers = [
        '[CH]=CC(=C)[CH]O[O](20187)',
    ],
    reactants = [
        ('C2H2(1342)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4979',
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


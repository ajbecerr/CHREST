species(
    label = '[CH2][C]=CC(=C)[CH]O[O](20233)',
    structure = SMILES('[CH2]C([CH]O[O])=C[C]=C'),
    E0 = (523.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,378.261,378.741],'cm^-1')),
        HinderedRotor(inertia=(0.658427,'amu*angstrom^2'), symmetry=1, barrier=(66.862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658247,'amu*angstrom^2'), symmetry=1, barrier=(66.8569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657354,'amu*angstrom^2'), symmetry=1, barrier=(66.8613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658181,'amu*angstrom^2'), symmetry=1, barrier=(66.8644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515395,0.0684685,-5.13794e-05,1.16843e-08,2.46331e-12,63102.7,30.0921], Tmin=(100,'K'), Tmax=(957.443,'K')), NASAPolynomial(coeffs=[16.0939,0.0218916,-7.40334e-06,1.25281e-09,-8.44813e-14,59271.3,-48.8178], Tmin=(957.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
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
    label = '[CH2]C([CH][C]=C)=C[O](25046)',
    structure = SMILES('[CH2][C]=CC([CH2])=C[O]'),
    E0 = (460.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67809,'amu*angstrom^2'), symmetry=1, barrier=(38.5826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6769,'amu*angstrom^2'), symmetry=1, barrier=(38.5552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6732,'amu*angstrom^2'), symmetry=1, barrier=(38.4703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.819844,0.0544109,-7.4317e-06,-4.27959e-08,2.44202e-11,55498.8,25.3985], Tmin=(100,'K'), Tmax=(923.114,'K')), NASAPolynomial(coeffs=[20.2652,0.0104351,-1.43317e-06,1.46184e-10,-1.23671e-14,50192.4,-76.1553], Tmin=(923.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[CH]=C([CH2])C=[C][CH2](16906)',
    structure = SMILES('[CH]C([CH2])=C[C]=C'),
    E0 = (706.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,218.576,218.577,218.577,218.577],'cm^-1')),
        HinderedRotor(inertia=(1.49931,'amu*angstrom^2'), symmetry=1, barrier=(50.831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49932,'amu*angstrom^2'), symmetry=1, barrier=(50.831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49932,'amu*angstrom^2'), symmetry=1, barrier=(50.831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.87,'J/mol'), sigma=(5.74688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.46 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39391,0.0500702,-1.91473e-05,-1.0955e-08,8.38317e-12,85081.9,22.1753], Tmin=(100,'K'), Tmax=(944.297,'K')), NASAPolynomial(coeffs=[11.1102,0.0263477,-9.15985e-06,1.54656e-09,-1.03096e-13,82469.5,-28.2591], Tmin=(944.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2][C]=C[C]=C(17191)',
    structure = SMILES('[CH2][C]=C[C]=C'),
    E0 = (612.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(2.17267,'amu*angstrom^2'), symmetry=1, barrier=(49.9539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16922,'amu*angstrom^2'), symmetry=1, barrier=(49.8747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20155,0.0347757,-1.38759e-05,-8.03191e-09,6.4918e-12,73770.6,17.1921], Tmin=(100,'K'), Tmax=(915.397,'K')), NASAPolynomial(coeffs=[9.16602,0.0169746,-5.40492e-06,8.73114e-10,-5.70733e-14,71966.3,-18.6825], Tmin=(915.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C)"""),
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
    label = 'C=[C]C=[C][CH]O[O](28485)',
    structure = SMILES('[CH2][C]=C[C]=CO[O]'),
    E0 = (680.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(2.5685,'amu*angstrom^2'), symmetry=1, barrier=(59.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767157,'amu*angstrom^2'), symmetry=1, barrier=(17.6385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.56822,'amu*angstrom^2'), symmetry=1, barrier=(59.0485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01922,0.0637119,-7.75705e-05,4.9893e-08,-1.25258e-11,81953.6,26.7867], Tmin=(100,'K'), Tmax=(1022.57,'K')), NASAPolynomial(coeffs=[12.5937,0.0164116,-5.21691e-06,7.8619e-10,-4.67953e-14,79692.3,-28.7944], Tmin=(1022.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
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
    label = '[CH]C(=C)[CH]O[O](22675)',
    structure = SMILES('[CH]C(=C)[CH]O[O]'),
    E0 = (525.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,426.391,427.662,428.94,429.241],'cm^-1')),
        HinderedRotor(inertia=(0.407888,'amu*angstrom^2'), symmetry=1, barrier=(51.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395576,'amu*angstrom^2'), symmetry=1, barrier=(51.0862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39968,'amu*angstrom^2'), symmetry=1, barrier=(50.9909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57252,0.0466541,-2.84193e-05,4.36067e-09,1.43993e-12,63286.8,23.6191], Tmin=(100,'K'), Tmax=(1077.55,'K')), NASAPolynomial(coeffs=[11.3392,0.0213978,-8.57272e-06,1.55485e-09,-1.06925e-13,60543.4,-27.1922], Tmin=(1077.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C([C]O[O])=C[C]=C(28486)',
    structure = SMILES('[CH2][C]=CC([CH2])=[C]O[O]'),
    E0 = (835.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,234.525],'cm^-1')),
        HinderedRotor(inertia=(0.20383,'amu*angstrom^2'), symmetry=1, barrier=(7.85371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00302488,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09078,'amu*angstrom^2'), symmetry=1, barrier=(82.6902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12419,'amu*angstrom^2'), symmetry=1, barrier=(82.4947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678048,0.0737998,-8.76725e-05,5.62684e-08,-1.44219e-11,100565,32.4313], Tmin=(100,'K'), Tmax=(953.215,'K')), NASAPolynomial(coeffs=[12.4896,0.0242348,-9.67626e-06,1.71905e-09,-1.15297e-13,98312.9,-23.9869], Tmin=(953.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]C([CH]O[O])=C[C]=C(20221)',
    structure = SMILES('[CH]C([CH]O[O])=C[C]=C'),
    E0 = (776.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192024,0.0777896,-7.58627e-05,3.89187e-08,-7.96183e-12,93496.2,30.8041], Tmin=(100,'K'), Tmax=(1186.92,'K')), NASAPolynomial(coeffs=[15.7788,0.0252606,-9.47717e-06,1.63099e-09,-1.07878e-13,89796.2,-47.0639], Tmin=(1186.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(776.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]=CC(=C)C1OO1(23665)',
    structure = SMILES('[CH2]C(=C[C]=C)C1OO1'),
    E0 = (370.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712001,0.0587941,-1.33806e-05,-3.78552e-08,2.34638e-11,44637.2,26.9283], Tmin=(100,'K'), Tmax=(898.032,'K')), NASAPolynomial(coeffs=[18.9646,0.0148185,-2.27183e-06,1.8045e-10,-9.16522e-15,39853.8,-67.5468], Tmin=(898.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(dioxirane) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=C1CC1O[O](28487)',
    structure = SMILES('C=C=C[C]1CC1O[O]'),
    E0 = (458.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07909,0.0515911,-9.3611e-06,-2.69948e-08,1.44441e-11,55262.9,27.0489], Tmin=(100,'K'), Tmax=(992.794,'K')), NASAPolynomial(coeffs=[15.5662,0.0220584,-8.30919e-06,1.55541e-09,-1.12482e-13,50965.2,-49.8957], Tmin=(992.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(ROOJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C1=CC(=C)C1O[O](22695)',
    structure = SMILES('C=C1[CH]C(=C)C1O[O]'),
    E0 = (343.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52685,0.0325152,5.70761e-05,-1.03131e-07,4.32144e-11,41427,26.1681], Tmin=(100,'K'), Tmax=(951.355,'K')), NASAPolynomial(coeffs=[18.331,0.0167967,-4.75668e-06,8.9576e-10,-7.20511e-14,35743.6,-67.1299], Tmin=(951.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C1=CC(=CO[O])C1(22672)',
    structure = SMILES('C=C1[CH]C(=CO[O])C1'),
    E0 = (365.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5569,0.0298805,6.82429e-05,-1.19897e-07,5.10099e-11,44077.6,26.8163], Tmin=(100,'K'), Tmax=(932.158,'K')), NASAPolynomial(coeffs=[19.8543,0.012329,-1.61672e-06,2.27942e-10,-2.399e-14,38017.7,-74.38], Tmin=(932.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2]C(=C=C=C)CO[O](28488)',
    structure = SMILES('[CH2]C#CC(=C)CO[O]'),
    E0 = (414.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12628,0.0670892,-7.11831e-05,4.48085e-08,-1.19838e-11,49995.2,28.5436], Tmin=(100,'K'), Tmax=(889.23,'K')), NASAPolynomial(coeffs=[8.4934,0.033951,-1.52856e-05,2.90285e-09,-2.02737e-13,48684.9,-6.13385], Tmin=(889.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(ROOJ) + radical(Propargyl)"""),
)

species(
    label = 'C=C=C=C(C)[CH]O[O](28489)',
    structure = SMILES('[CH2]C#CC(C)=CO[O]'),
    E0 = (412.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910315,0.0681605,-6.87632e-05,3.79995e-08,-8.57884e-12,49685.4,28.6945], Tmin=(100,'K'), Tmax=(1063.62,'K')), NASAPolynomial(coeffs=[11.7015,0.0275768,-1.15278e-05,2.12424e-09,-1.46331e-13,47389.9,-24.0323], Tmin=(1063.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(ROOJ) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C]1C([C]=C)C1O[O](28490)',
    structure = SMILES('[CH2][C]1C([C]=C)C1O[O]'),
    E0 = (766.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32534,0.0486057,-1.06563e-05,-2.22552e-08,1.25892e-11,92279.1,33.4811], Tmin=(100,'K'), Tmax=(968.818,'K')), NASAPolynomial(coeffs=[13.3387,0.0229427,-7.98443e-06,1.40916e-09,-9.82894e-14,88827.9,-29.8937], Tmin=(968.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(766.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C[C]1[CH]O[O](28491)',
    structure = SMILES('C=[C]C1C[C]1[CH]O[O]'),
    E0 = (762.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23402,0.0526791,-2.83701e-05,2.28521e-09,2.02795e-12,91816.1,33.2205], Tmin=(100,'K'), Tmax=(1129.49,'K')), NASAPolynomial(coeffs=[12.2424,0.02656,-1.077e-05,1.98233e-09,-1.37283e-13,88508.6,-24.8623], Tmin=(1129.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC1([CH2])[CH]O[O](22938)',
    structure = SMILES('[CH2]C1([CH]O[O])[CH]C1=C'),
    E0 = (720.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.591614,0.0688089,-5.95203e-05,2.6047e-08,-4.56151e-12,86827.4,26.8816], Tmin=(100,'K'), Tmax=(1364,'K')), NASAPolynomial(coeffs=[15.9365,0.0238089,-1.00332e-05,1.85955e-09,-1.28304e-13,82641.3,-51.9118], Tmin=(1364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(720.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Allyl_S) + radical(CCsJOOH) + radical(Neopentyl)"""),
)

species(
    label = '[CH2][C]1[CH]OOC1[C]=C(28492)',
    structure = SMILES('[CH2][C]1[CH]OOC1[C]=C'),
    E0 = (730.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59163,0.0379974,3.26461e-05,-7.88888e-08,3.72385e-11,87932.8,31.7944], Tmin=(100,'K'), Tmax=(884.09,'K')), NASAPolynomial(coeffs=[15.2782,0.0172862,-2.13812e-06,6.8689e-11,9.82729e-16,83902.2,-41.6583], Tmin=(884.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(730.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(C2CJCOOH) + radical(CCsJOOC) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1([CH][C]=C)[CH]OO1(28493)',
    structure = SMILES('[CH2]C1([CH][C]=C)[CH]OO1'),
    E0 = (742.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.427987,0.0821241,-9.3955e-05,5.82764e-08,-1.47352e-11,89426.8,24.6796], Tmin=(100,'K'), Tmax=(954.053,'K')), NASAPolynomial(coeffs=[12.4729,0.0316236,-1.45557e-05,2.79405e-09,-1.96509e-13,87128.5,-32.8638], Tmin=(954.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(742.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=CCJCO) + radical(CCsJOO) + radical(CJCOOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C=C=C)[CH]O[O](28494)',
    structure = SMILES('[CH2]C#CC(=C)[CH]O[O]'),
    E0 = (532.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2100,2250,500,550,847.053,847.902],'cm^-1')),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235142,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553354,'amu*angstrom^2'), symmetry=1, barrier=(12.7227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538818,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.795448,0.0624418,-5.37629e-05,2.3481e-08,-4.07631e-12,64124.5,30.3223], Tmin=(100,'K'), Tmax=(1386.37,'K')), NASAPolynomial(coeffs=[15.5045,0.0200028,-7.84546e-06,1.40058e-09,-9.46242e-14,60046,-45.4457], Tmin=(1386.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(ROOJ) + radical(C=CCJO) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=CC(=C)[CH]O[O](22639)',
    structure = SMILES('C#CC=C([CH2])[CH]O[O]'),
    E0 = (504.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,400.959],'cm^-1')),
        HinderedRotor(inertia=(0.29897,'amu*angstrom^2'), symmetry=1, barrier=(34.1077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677127,'amu*angstrom^2'), symmetry=1, barrier=(77.2496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230901,'amu*angstrom^2'), symmetry=1, barrier=(26.3422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677128,'amu*angstrom^2'), symmetry=1, barrier=(77.2496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46623,0.0700507,-6.80174e-05,3.35713e-08,-6.51441e-12,60836.6,27.6612], Tmin=(100,'K'), Tmax=(1257.53,'K')), NASAPolynomial(coeffs=[16.555,0.018875,-6.97419e-06,1.20985e-09,-8.08567e-14,56790.1,-53.6443], Tmin=(1257.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(ROOJ) + radical(C=CCJO) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH2]C(C=O)=C[C]=C(28495)',
    structure = SMILES('[CH2]C(C=C=C)=C[O]'),
    E0 = (281.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43467,'amu*angstrom^2'), symmetry=1, barrier=(32.9858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4326,'amu*angstrom^2'), symmetry=1, barrier=(32.9383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612329,0.0583586,-1.49765e-05,-3.75896e-08,2.2973e-11,33942.8,23.1651], Tmin=(100,'K'), Tmax=(929.866,'K')), NASAPolynomial(coeffs=[21.8818,0.00821495,-7.93005e-07,6.55667e-11,-8.55915e-15,28199.5,-87.5141], Tmin=(929.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C=C=CO[O](28496)',
    structure = SMILES('C=[C]C=C=CO[O]'),
    E0 = (501.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.17931,'amu*angstrom^2'), symmetry=1, barrier=(27.1147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17738,'amu*angstrom^2'), symmetry=1, barrier=(27.0703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.920974,0.0663002,-8.00166e-05,4.79639e-08,-1.0681e-11,60393,24.1659], Tmin=(100,'K'), Tmax=(879.114,'K')), NASAPolynomial(coeffs=[13.9309,0.0146798,-4.86382e-06,7.7433e-10,-4.87492e-14,57812.8,-38.5882], Tmin=(879.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=C=C(17211)',
    structure = SMILES('C=[C]C=C=C'),
    E0 = (433.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.5987,'amu*angstrom^2'), symmetry=1, barrier=(36.7571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98554,0.0388285,-2.1811e-05,-2.2868e-09,4.80017e-12,52215,14.9887], Tmin=(100,'K'), Tmax=(940.235,'K')), NASAPolynomial(coeffs=[10.8095,0.0147077,-4.73738e-06,7.85956e-10,-5.27182e-14,49962.5,-30.1922], Tmin=(940.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(=[C][C]=C)[CH]O[O](28497)',
    structure = SMILES('[CH2]C(=[C][C]=C)[CH]O[O]'),
    E0 = (722.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,324.871,326.266],'cm^-1')),
        HinderedRotor(inertia=(0.75093,'amu*angstrom^2'), symmetry=1, barrier=(56.3966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.738139,'amu*angstrom^2'), symmetry=1, barrier=(56.3188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743711,'amu*angstrom^2'), symmetry=1, barrier=(56.2359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.736918,'amu*angstrom^2'), symmetry=1, barrier=(56.2401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482908,0.0744428,-7.83255e-05,4.15402e-08,-7.90446e-12,87032.8,29.8077], Tmin=(100,'K'), Tmax=(899.109,'K')), NASAPolynomial(coeffs=[14.4437,0.0221721,-7.536e-06,1.22252e-09,-7.8167e-14,84124.7,-38.2721], Tmin=(899.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C][CH]C(=C)[CH]O[O](22641)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH]O[O]'),
    E0 = (770.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,314.252],'cm^-1')),
        HinderedRotor(inertia=(0.890749,'amu*angstrom^2'), symmetry=1, barrier=(62.4201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.89073,'amu*angstrom^2'), symmetry=1, barrier=(62.4201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890704,'amu*angstrom^2'), symmetry=1, barrier=(62.4199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890722,'amu*angstrom^2'), symmetry=1, barrier=(62.4201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446445,0.0719331,-6.78292e-05,2.87929e-08,-3.33456e-12,92822.2,30.7429], Tmin=(100,'K'), Tmax=(951.677,'K')), NASAPolynomial(coeffs=[16.3184,0.0190858,-6.38491e-06,1.05773e-09,-6.992e-14,89173.4,-48.3423], Tmin=(951.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(770.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=[C]C=C)[CH]O[O](28498)',
    structure = SMILES('[CH2]C(=[C]C=C)[CH]O[O]'),
    E0 = (523.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,378.15,378.195],'cm^-1')),
        HinderedRotor(inertia=(0.657393,'amu*angstrom^2'), symmetry=1, barrier=(66.8608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658514,'amu*angstrom^2'), symmetry=1, barrier=(66.865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.657618,'amu*angstrom^2'), symmetry=1, barrier=(66.8637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.656345,'amu*angstrom^2'), symmetry=1, barrier=(66.8518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515395,0.0684685,-5.13794e-05,1.16843e-08,2.46331e-12,63102.7,30.0921], Tmin=(100,'K'), Tmax=(957.443,'K')), NASAPolynomial(coeffs=[16.0939,0.0218916,-7.40334e-06,1.25281e-09,-8.44813e-14,59271.3,-48.8178], Tmin=(957.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[CH]C(=C)[CH]O[O](22642)',
    structure = SMILES('[CH]=CC=C([CH2])[CH]O[O]'),
    E0 = (571.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,414.213],'cm^-1')),
        HinderedRotor(inertia=(0.200844,'amu*angstrom^2'), symmetry=1, barrier=(24.4295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631605,'amu*angstrom^2'), symmetry=1, barrier=(76.2644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200955,'amu*angstrom^2'), symmetry=1, barrier=(24.4274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62787,'amu*angstrom^2'), symmetry=1, barrier=(76.2456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521768,0.0654208,-3.88406e-05,-3.9578e-09,8.38535e-12,68890.2,30.8758], Tmin=(100,'K'), Tmax=(963.466,'K')), NASAPolynomial(coeffs=[17.8649,0.0189882,-6.36045e-06,1.11406e-09,-7.8424e-14,64361.5,-58.3083], Tmin=(963.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=[C]C(=C)CO[O](20234)',
    structure = SMILES('[CH2]C(=[C][C]=C)CO[O]'),
    E0 = (605.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,562.897],'cm^-1')),
        HinderedRotor(inertia=(4.07605,'amu*angstrom^2'), symmetry=1, barrier=(93.7164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.07418,'amu*angstrom^2'), symmetry=1, barrier=(93.6734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.07742,'amu*angstrom^2'), symmetry=1, barrier=(93.748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.07674,'amu*angstrom^2'), symmetry=1, barrier=(93.7324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331225,0.0848602,-0.000116214,8.96598e-08,-2.73011e-11,72924.2,29.7536], Tmin=(100,'K'), Tmax=(913.329,'K')), NASAPolynomial(coeffs=[10.3408,0.0312837,-1.22287e-05,2.08265e-09,-1.3343e-13,71502,-15.4055], Tmin=(913.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C][C]=C(C)[CH]O[O](28499)',
    structure = SMILES('C=[C][C]=C(C)[CH]O[O]'),
    E0 = (604.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.775242,'amu*angstrom^2'), symmetry=1, barrier=(17.8243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82322,'amu*angstrom^2'), symmetry=1, barrier=(64.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82106,'amu*angstrom^2'), symmetry=1, barrier=(64.8616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82179,'amu*angstrom^2'), symmetry=1, barrier=(64.8785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259002,0.0826715,-9.92031e-05,6.44482e-08,-1.66722e-11,72839,29.0547], Tmin=(100,'K'), Tmax=(946.481,'K')), NASAPolynomial(coeffs=[13.5441,0.0265255,-1.02207e-05,1.77129e-09,-1.16753e-13,70324.2,-34.3072], Tmin=(946.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(604.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C][CH]C(=C)CO[O](20241)',
    structure = SMILES('[CH]=[C]C=C([CH2])CO[O]'),
    E0 = (653.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(3.63112,'amu*angstrom^2'), symmetry=1, barrier=(83.4865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.629,'amu*angstrom^2'), symmetry=1, barrier=(83.4379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352395,'amu*angstrom^2'), symmetry=1, barrier=(8.10225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.62962,'amu*angstrom^2'), symmetry=1, barrier=(83.452,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.487694,0.0798819,-9.61005e-05,6.28509e-08,-1.59225e-11,78705.4,30.0088], Tmin=(100,'K'), Tmax=(775.404,'K')), NASAPolynomial(coeffs=[11.7158,0.0290811,-1.16018e-05,2.04425e-09,-1.35828e-13,76750.1,-22.6843], Tmin=(775.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(C)=CO[O](22653)',
    structure = SMILES('[CH]=C=C[C](C)[CH]O[O]'),
    E0 = (651.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,540,610,2055,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.16138,'amu*angstrom^2'), symmetry=1, barrier=(3.71045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160795,'amu*angstrom^2'), symmetry=1, barrier=(3.697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.88435,'amu*angstrom^2'), symmetry=1, barrier=(66.3168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87949,'amu*angstrom^2'), symmetry=1, barrier=(66.2051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666102,0.0763,-9.05373e-05,6.14582e-08,-1.70849e-11,78489.2,30.1063], Tmin=(100,'K'), Tmax=(872.114,'K')), NASAPolynomial(coeffs=[10.4526,0.0314126,-1.33312e-05,2.43854e-09,-1.66009e-13,76782.2,-15.7687], Tmin=(872.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCJ(C)CO) + radical(CCsJOOH) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2]C(=[C][C]=C)[CH]OO(28500)',
    structure = SMILES('[CH2]C(=[C][C]=C)[CH]OO'),
    E0 = (570.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,305.994,306.731],'cm^-1')),
        HinderedRotor(inertia=(0.768027,'amu*angstrom^2'), symmetry=1, barrier=(50.9504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.761492,'amu*angstrom^2'), symmetry=1, barrier=(50.9784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773835,'amu*angstrom^2'), symmetry=1, barrier=(50.9768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.761402,'amu*angstrom^2'), symmetry=1, barrier=(50.9252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773014,'amu*angstrom^2'), symmetry=1, barrier=(50.977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125126,0.0791747,-7.77568e-05,3.78013e-08,-6.61264e-12,68766.4,30.3304], Tmin=(100,'K'), Tmax=(984.109,'K')), NASAPolynomial(coeffs=[16.5103,0.0230075,-8.04625e-06,1.34879e-09,-8.87452e-14,65036.3,-51.0229], Tmin=(984.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C][CH]C(=C)[CH]OO(22643)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH]OO'),
    E0 = (618.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,311.59],'cm^-1')),
        HinderedRotor(inertia=(0.675409,'amu*angstrom^2'), symmetry=1, barrier=(46.533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.67541,'amu*angstrom^2'), symmetry=1, barrier=(46.533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675411,'amu*angstrom^2'), symmetry=1, barrier=(46.533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675412,'amu*angstrom^2'), symmetry=1, barrier=(46.533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.67541,'amu*angstrom^2'), symmetry=1, barrier=(46.533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.143002,0.0759871,-6.47094e-05,2.14744e-08,-3.88103e-13,74553.4,31.073], Tmin=(100,'K'), Tmax=(977.763,'K')), NASAPolynomial(coeffs=[18.2412,0.0201728,-7.0433e-06,1.21952e-09,-8.34772e-14,70143.1,-60.2885], Tmin=(977.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CC([C]=C)O[O](21296)',
    structure = SMILES('C=[C][CH]C([C]=C)O[O]'),
    E0 = (660.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,376.897,376.905,376.921],'cm^-1')),
        HinderedRotor(inertia=(0.00118675,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27058,'amu*angstrom^2'), symmetry=1, barrier=(27.2787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270611,'amu*angstrom^2'), symmetry=1, barrier=(27.2793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270563,'amu*angstrom^2'), symmetry=1, barrier=(27.2789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3940.4,'J/mol'), sigma=(6.57992,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.48 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366681,0.06927,-5.96723e-05,2.54595e-08,-4.27754e-12,79547,34.6298], Tmin=(100,'K'), Tmax=(1437.8,'K')), NASAPolynomial(coeffs=[18.2373,0.0195533,-7.80471e-06,1.41e-09,-9.58787e-14,74408.1,-58.0749], Tmin=(1437.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])[CH]O[O](28501)',
    structure = SMILES('[CH2][C]=[C]C([CH2])[CH]O[O]'),
    E0 = (1000.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.157053,'amu*angstrom^2'), symmetry=1, barrier=(3.61095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157041,'amu*angstrom^2'), symmetry=1, barrier=(3.61067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157083,'amu*angstrom^2'), symmetry=1, barrier=(3.61166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.76848,'amu*angstrom^2'), symmetry=1, barrier=(63.6527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023951,'amu*angstrom^2'), symmetry=1, barrier=(63.6503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.180497,0.0921413,-0.000143991,1.22578e-07,-4.04655e-11,120430,34.9763], Tmin=(100,'K'), Tmax=(881.646,'K')), NASAPolynomial(coeffs=[9.33735,0.0333721,-1.46978e-05,2.65145e-09,-1.75296e-13,119485,-4.25017], Tmin=(881.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1000.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]1CC1O[O](28502)',
    structure = SMILES('[CH2][C]=C[C]1CC1O[O]'),
    E0 = (671.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00222,0.0541035,-1.67211e-05,-1.93068e-08,1.17299e-11,80851.2,29.0254], Tmin=(100,'K'), Tmax=(998.051,'K')), NASAPolynomial(coeffs=[15.4944,0.0223156,-8.46462e-06,1.57525e-09,-1.1301e-13,76648.8,-47.4237], Tmin=(998.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=C[C]([CH]O[O])C1(28503)',
    structure = SMILES('[CH2][C]1C=C([CH]O[O])C1'),
    E0 = (598.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32721,0.0432492,1.89979e-05,-6.0776e-08,2.81267e-11,72080.4,28.251], Tmin=(100,'K'), Tmax=(939.649,'K')), NASAPolynomial(coeffs=[16.0569,0.0196142,-5.63849e-06,9.50981e-10,-6.86323e-14,67587.6,-51.0721], Tmin=(939.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(Allyl_T) + radical(C=CCJO) + radical(Isobutyl)"""),
)

species(
    label = '[O]O[CH][C]1C=[C]CC1(28504)',
    structure = SMILES('[O]O[CH]C1=C[C]CC1'),
    E0 = (618.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20897,0.0442699,1.9283e-05,-6.0376e-08,2.69397e-11,74543.5,27.218], Tmin=(100,'K'), Tmax=(974.273,'K')), NASAPolynomial(coeffs=[17.1939,0.0200586,-7.20689e-06,1.3834e-09,-1.04231e-13,69463.1,-59.5713], Tmin=(974.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(ROOJ) + radical(C=CCJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[C]1[CH]OOC1(28505)',
    structure = SMILES('[CH2][C]=C[C]1[CH]OOC1'),
    E0 = (640.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28629,0.0455067,1.67e-05,-6.29434e-08,3.11461e-11,77104.5,26.4633], Tmin=(100,'K'), Tmax=(894.727,'K')), NASAPolynomial(coeffs=[15.8486,0.0190046,-3.58274e-06,3.87683e-10,-2.24333e-14,72953.5,-50.8056], Tmin=(894.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(12dioxolane) + radical(CCJ(C)CO) + radical(CCsJOOC) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1([CH]O[O])C=[C]C1(22637)',
    structure = SMILES('[CH2]C1([CH]O[O])C=[C]C1'),
    E0 = (789.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89073,0.0656584,-5.56165e-05,2.42047e-08,-4.28759e-12,95100,28.5086], Tmin=(100,'K'), Tmax=(1329.79,'K')), NASAPolynomial(coeffs=[13.8142,0.0267846,-1.17668e-05,2.22138e-09,-1.54717e-13,91662.9,-37.523], Tmin=(1329.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(789.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(CCsJOOH) + radical(Neopentyl) + radical(cyclobutene-vinyl)"""),
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
    label = '[CH]=C(C=[C][CH2])CO[O](19510)',
    structure = SMILES('[CH]C(=C[C]=C)CO[O]'),
    E0 = (658.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4016.09,'J/mol'), sigma=(6.63908,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=627.30 K, Pc=31.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273025,0.0853672,-0.000103427,7.32099e-08,-2.13125e-11,79377.8,29.9223], Tmin=(100,'K'), Tmax=(833.417,'K')), NASAPolynomial(coeffs=[10.471,0.0364288,-1.53598e-05,2.77297e-09,-1.8657e-13,77677.7,-17.4203], Tmin=(833.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(658.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([C]=[C]C)[CH]O[O](28506)',
    structure = SMILES('C=C([C]=[C]C)[CH]O[O]'),
    E0 = (644.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,202.086,202.088],'cm^-1')),
        HinderedRotor(inertia=(0.536981,'amu*angstrom^2'), symmetry=1, barrier=(15.5612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.02549,'amu*angstrom^2'), symmetry=1, barrier=(87.6775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.02535,'amu*angstrom^2'), symmetry=1, barrier=(87.6774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536986,'amu*angstrom^2'), symmetry=1, barrier=(15.5613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.391379,0.0796758,-9.1123e-05,5.58168e-08,-1.3672e-11,77682.6,29.3023], Tmin=(100,'K'), Tmax=(994.816,'K')), NASAPolynomial(coeffs=[13.6542,0.0263489,-1.07168e-05,1.93404e-09,-1.31348e-13,75043.8,-34.6144], Tmin=(994.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]C([CH]O[O])=CC=C(20227)',
    structure = SMILES('[CH]=C([CH]O[O])C=C[CH2]'),
    E0 = (573.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,209.769],'cm^-1')),
        HinderedRotor(inertia=(0.00394002,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33573,'amu*angstrom^2'), symmetry=1, barrier=(41.8829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36357,'amu*angstrom^2'), symmetry=1, barrier=(41.852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40092,'amu*angstrom^2'), symmetry=1, barrier=(41.935,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521768,0.0654208,-3.88406e-05,-3.95778e-09,8.38534e-12,69066.3,30.7198], Tmin=(100,'K'), Tmax=(963.466,'K')), NASAPolynomial(coeffs=[17.8649,0.0189882,-6.36045e-06,1.11406e-09,-7.8424e-14,64537.6,-58.4643], Tmin=(963.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C=[C]C(20236)',
    structure = SMILES('[CH]=C([CH]O[O])C=[C]C'),
    E0 = (692.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,220.715],'cm^-1')),
        HinderedRotor(inertia=(0.454104,'amu*angstrom^2'), symmetry=1, barrier=(15.6931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453459,'amu*angstrom^2'), symmetry=1, barrier=(15.684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.789811,'amu*angstrom^2'), symmetry=1, barrier=(27.3673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454385,'amu*angstrom^2'), symmetry=1, barrier=(15.6877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336182,0.0773977,-8.14361e-05,4.40166e-08,-9.41844e-12,83472.7,30.3034], Tmin=(100,'K'), Tmax=(1137.78,'K')), NASAPolynomial(coeffs=[15.9271,0.0225864,-9.17597e-06,1.67716e-09,-1.1546e-13,79924.9,-46.9266], Tmin=(1137.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(692.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C([CH]OO)=C[C]=C(20226)',
    structure = SMILES('[CH]C([CH]OO)=C[C]=C'),
    E0 = (624.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.228982,0.0831545,-7.69579e-05,3.65119e-08,-6.86341e-12,75232.6,31.5609], Tmin=(100,'K'), Tmax=(1291.14,'K')), NASAPolynomial(coeffs=[18.5463,0.0249882,-9.38247e-06,1.62013e-09,-1.07438e-13,70384.3,-63.8167], Tmin=(1291.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CCJO) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C1C(=C)C1O[O](22556)',
    structure = SMILES('C=[C]C1C(=C)C1O[O]'),
    E0 = (547.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787709,0.0693833,-6.63414e-05,3.34422e-08,-6.83525e-12,65928.6,25.2222], Tmin=(100,'K'), Tmax=(1170.65,'K')), NASAPolynomial(coeffs=[13.3463,0.0264721,-1.1358e-05,2.13012e-09,-1.48385e-13,62988.2,-37.3443], Tmin=(1170.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([CH]O[O])C[C]=C(20080)',
    structure = SMILES('[CH]=C([CH]O[O])C[C]=C'),
    E0 = (722.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,563.943],'cm^-1')),
        HinderedRotor(inertia=(0.010834,'amu*angstrom^2'), symmetry=1, barrier=(2.44035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064499,'amu*angstrom^2'), symmetry=1, barrier=(14.5172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50385,'amu*angstrom^2'), symmetry=1, barrier=(34.5765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631387,'amu*angstrom^2'), symmetry=1, barrier=(14.5168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.644475,0.0716765,-6.98211e-05,3.5401e-08,-7.22043e-12,86995.4,32.3759], Tmin=(100,'K'), Tmax=(1178.73,'K')), NASAPolynomial(coeffs=[14.3121,0.0252958,-1.07995e-05,2.01971e-09,-1.40535e-13,83773.3,-35.81], Tmin=(1178.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(=C)[CH]O[O](20307)',
    structure = SMILES('[CH]=[C]CC(=C)[CH]O[O]'),
    E0 = (722.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,563.943],'cm^-1')),
        HinderedRotor(inertia=(0.010834,'amu*angstrom^2'), symmetry=1, barrier=(2.44035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.064499,'amu*angstrom^2'), symmetry=1, barrier=(14.5172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50385,'amu*angstrom^2'), symmetry=1, barrier=(34.5765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631387,'amu*angstrom^2'), symmetry=1, barrier=(14.5168,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3938.88,'J/mol'), sigma=(6.5717,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.24 K, Pc=31.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.644475,0.0716765,-6.98211e-05,3.5401e-08,-7.22043e-12,86995.4,32.3759], Tmin=(100,'K'), Tmax=(1178.73,'K')), NASAPolynomial(coeffs=[14.3121,0.0252958,-1.07995e-05,2.01971e-09,-1.40535e-13,83773.3,-35.81], Tmin=(1178.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C=CC[C][CH]O[O](20294)',
    structure = SMILES('C=C=CC[C][CH]O[O]'),
    E0 = (800.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,610,2055,239.46,240.144,241.858,242.671],'cm^-1')),
        HinderedRotor(inertia=(1.82485,'amu*angstrom^2'), symmetry=1, barrier=(75.8944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188554,'amu*angstrom^2'), symmetry=1, barrier=(7.8712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694696,'amu*angstrom^2'), symmetry=1, barrier=(28.4667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191523,'amu*angstrom^2'), symmetry=1, barrier=(7.87235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4003.42,'J/mol'), sigma=(6.71817,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=625.32 K, Pc=29.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25579,0.0886339,-0.000128569,1.05368e-07,-3.48359e-11,96437.6,30.902], Tmin=(100,'K'), Tmax=(803.661,'K')), NASAPolynomial(coeffs=[10.0366,0.0341497,-1.60459e-05,3.04157e-09,-2.09491e-13,95052.9,-12.981], Tmin=(803.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]=CC1=COOC1(19223)',
    structure = SMILES('C=[C]C=C1[CH]OOC1'),
    E0 = (346.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44536,0.0336811,5.93083e-05,-1.07195e-07,4.50265e-11,41791.6,25.5394], Tmin=(100,'K'), Tmax=(948.254,'K')), NASAPolynomial(coeffs=[18.6672,0.0176497,-4.8891e-06,9.0083e-10,-7.18896e-14,35980.1,-70.0523], Tmin=(948.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1=COOC([CH2])=C1(28507)',
    structure = SMILES('[CH2]C1=CC(=C)OO[CH]1'),
    E0 = (231.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6611,0.0261357,8.41867e-05,-1.26569e-07,4.8852e-11,27940.1,20.5741], Tmin=(100,'K'), Tmax=(987.305,'K')), NASAPolynomial(coeffs=[17.5246,0.0251082,-1.03356e-05,2.13554e-09,-1.66165e-13,21725.3,-71.3652], Tmin=(987.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJO) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=CC1=COC1(20232)',
    structure = SMILES('C=[C]C=C1[CH]OC1'),
    E0 = (371.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81841,0.0374024,9.85119e-06,-3.56906e-08,1.52513e-11,44740.9,22.6671], Tmin=(100,'K'), Tmax=(1005.35,'K')), NASAPolynomial(coeffs=[10.7885,0.026268,-1.01727e-05,1.88222e-09,-1.33274e-13,41696.4,-26.8278], Tmin=(1005.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJ(O)C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1=COC([CH2])=C1(28508)',
    structure = SMILES('[CH2]C1=CC(=C)O[CH]1'),
    E0 = (145.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1025,0.015718,9.95797e-05,-1.49159e-07,6.09031e-11,17621.1,19.5871], Tmin=(100,'K'), Tmax=(929.263,'K')), NASAPolynomial(coeffs=[18.6237,0.0112647,-8.36237e-07,7.75764e-11,-1.45348e-14,11672.3,-74.3925], Tmin=(929.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJ(O)C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]1C=C([CH2])C1O[O](28509)',
    structure = SMILES('[CH2][C]1C=C([CH2])C1O[O]'),
    E0 = (633.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04972,0.0495869,4.7105e-06,-4.8116e-08,2.38754e-11,76314.9,29.6908], Tmin=(100,'K'), Tmax=(949.736,'K')), NASAPolynomial(coeffs=[17.4439,0.0183158,-5.56391e-06,9.76967e-10,-7.1711e-14,71497.1,-57.5258], Tmin=(949.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]1C=[C]CC1O[O](28510)',
    structure = SMILES('C=C1[CH][C]CC1O[O]'),
    E0 = (633.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3606,0.0400812,3.01509e-05,-6.92241e-08,2.92087e-11,76263.4,24.7372], Tmin=(100,'K'), Tmax=(985.771,'K')), NASAPolynomial(coeffs=[16.5559,0.0219372,-8.45445e-06,1.66424e-09,-1.25741e-13,71153.3,-59.0776], Tmin=(985.771,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclopentane) + radical(ROOJ) + radical(Allyl_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[C]([CH2])C1OO1(28511)',
    structure = SMILES('[CH2][C]=C[C]([CH2])C1OO1'),
    E0 = (696.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898083,0.0541964,-5.33251e-06,-4.27892e-08,2.41049e-11,83936.5,31.8541], Tmin=(100,'K'), Tmax=(911.611,'K')), NASAPolynomial(coeffs=[18.1865,0.0155828,-3.08085e-06,3.81998e-10,-2.52323e-14,79236.8,-58.4408], Tmin=(911.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC([CH2])=[C]OO(28512)',
    structure = SMILES('[CH2][C]=CC([CH2])=[C]OO'),
    E0 = (683.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.0294273,'amu*angstrom^2'), symmetry=1, barrier=(14.9322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.39376,'amu*angstrom^2'), symmetry=1, barrier=(78.0291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648712,'amu*angstrom^2'), symmetry=1, barrier=(14.9152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154623,'amu*angstrom^2'), symmetry=1, barrier=(78.0595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.56929,'amu*angstrom^2'), symmetry=1, barrier=(78.024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419803,0.0773397,-8.28247e-05,4.67646e-08,-1.05485e-11,82293.9,32.5979], Tmin=(100,'K'), Tmax=(1077.05,'K')), NASAPolynomial(coeffs=[14.4426,0.025261,-1.02952e-05,1.87073e-09,-1.27964e-13,79273.2,-36.095], Tmin=(1077.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=CC(C)=[C]O[O](28513)',
    structure = SMILES('[CH2][C]=CC(C)=[C]O[O]'),
    E0 = (683.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678933,0.0762876,-9.6248e-05,6.98164e-08,-2.06252e-11,82341.1,32.148], Tmin=(100,'K'), Tmax=(824.672,'K')), NASAPolynomial(coeffs=[10.138,0.0303968,-1.27585e-05,2.30831e-09,-1.55504e-13,80781.3,-11.6611], Tmin=(824.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=CC([CH2])=[C]O[O](28514)',
    structure = SMILES('[CH2]C=CC([CH2])=[C]O[O]'),
    E0 = (597.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,267.545],'cm^-1')),
        HinderedRotor(inertia=(1.60742,'amu*angstrom^2'), symmetry=1, barrier=(80.0211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00235476,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23062,'amu*angstrom^2'), symmetry=1, barrier=(11.451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61673,'amu*angstrom^2'), symmetry=1, barrier=(79.9973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542916,0.071238,-7.0221e-05,3.6714e-08,-7.65164e-12,71968.5,32.4381], Tmin=(100,'K'), Tmax=(1166.26,'K')), NASAPolynomial(coeffs=[14.5172,0.0233089,-8.57575e-06,1.47553e-09,-9.78266e-14,68709,-37.1291], Tmin=(1166.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C=[C]C(28515)',
    structure = SMILES('[CH2]C(=[C]O[O])C=[C]C'),
    E0 = (717.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.47557,'amu*angstrom^2'), symmetry=1, barrier=(10.9343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.475576,'amu*angstrom^2'), symmetry=1, barrier=(10.9344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474964,'amu*angstrom^2'), symmetry=1, barrier=(10.9204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.475224,'amu*angstrom^2'), symmetry=1, barrier=(10.9263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.432161,0.0822788,-0.000109395,8.02656e-08,-2.36717e-11,86371.8,31.7579], Tmin=(100,'K'), Tmax=(829.01,'K')), NASAPolynomial(coeffs=[11.5404,0.0286823,-1.24196e-05,2.28246e-09,-1.5515e-13,84530,-19.75], Tmin=(829.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C1CC1=CO[O](28516)',
    structure = SMILES('C=[C]C1CC1=CO[O]'),
    E0 = (564.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723925,0.0702632,-7.06523e-05,3.78762e-08,-8.17141e-12,67976.5,26.7099], Tmin=(100,'K'), Tmax=(1119.77,'K')), NASAPolynomial(coeffs=[13.372,0.0250815,-1.01275e-05,1.84145e-09,-1.26132e-13,65143.9,-35.7403], Tmin=(1119.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=COOC1[C]=C(22696)',
    structure = SMILES('C=[C]C1OO[CH]C1=C'),
    E0 = (391.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59438,0.0290878,6.93106e-05,-1.13871e-07,4.579e-11,47199.7,27.5344], Tmin=(100,'K'), Tmax=(970.958,'K')), NASAPolynomial(coeffs=[18.4735,0.0189534,-6.80034e-06,1.39546e-09,-1.11737e-13,41121.9,-67.8191], Tmin=(970.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=COC1[C]=C(22697)',
    structure = SMILES('[CH2]C1=COC1[C]=C'),
    E0 = (413.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05154,0.0408098,4.10408e-05,-9.85983e-08,4.49657e-11,49804.8,22.4267], Tmin=(100,'K'), Tmax=(935.775,'K')), NASAPolynomial(coeffs=[23.8581,0.00453789,1.05796e-06,-2.07294e-10,3.62898e-15,42856.2,-100.409], Tmin=(935.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(=C)[C]O[O](22313)',
    structure = SMILES('[CH2]C(=[C]O[O])C[C]=C'),
    E0 = (746.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,299.38,306.569],'cm^-1')),
        HinderedRotor(inertia=(0.0018432,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125417,'amu*angstrom^2'), symmetry=1, barrier=(8.3207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123981,'amu*angstrom^2'), symmetry=1, barrier=(8.305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12632,'amu*angstrom^2'), symmetry=1, barrier=(8.26708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700687,0.0770123,-9.93452e-05,7.37434e-08,-2.24396e-11,89896.3,33.9744], Tmin=(100,'K'), Tmax=(798.494,'K')), NASAPolynomial(coeffs=[9.74283,0.0317096,-1.423e-05,2.66987e-09,-1.83992e-13,88452.5,-7.61267], Tmin=(798.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(746.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C][CH]O[O](22548)',
    structure = SMILES('[CH2][C][CH]O[O]'),
    E0 = (793.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1493.49],'cm^-1')),
        HinderedRotor(inertia=(0.167737,'amu*angstrom^2'), symmetry=1, barrier=(3.85659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168521,'amu*angstrom^2'), symmetry=1, barrier=(3.87462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16787,'amu*angstrom^2'), symmetry=1, barrier=(3.85967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07407,0.0474498,-7.98325e-05,7.1518e-08,-2.45548e-11,95458.7,21.5901], Tmin=(100,'K'), Tmax=(855.952,'K')), NASAPolynomial(coeffs=[6.67897,0.0159696,-7.20999e-06,1.35978e-09,-9.27035e-14,95035.3,2.22186], Tmin=(855.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(793.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=[C][C]([CH2])[CH]O[O](28517)',
    structure = SMILES('[CH2][CH][C]=C([CH2])[CH]O[O]'),
    E0 = (831.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,285.741,1224.73],'cm^-1')),
        HinderedRotor(inertia=(0.53165,'amu*angstrom^2'), symmetry=1, barrier=(30.7765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186462,'amu*angstrom^2'), symmetry=1, barrier=(19.8464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00449127,'amu*angstrom^2'), symmetry=1, barrier=(4.78036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0800366,'amu*angstrom^2'), symmetry=1, barrier=(85.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47276,'amu*angstrom^2'), symmetry=1, barrier=(85.1907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620422,0.0701322,-6.37036e-05,2.97931e-08,-5.60335e-12,100164,32.773], Tmin=(100,'K'), Tmax=(1272.52,'K')), NASAPolynomial(coeffs=[14.9568,0.0250671,-1.05819e-05,1.96262e-09,-1.35696e-13,96514.9,-39.8467], Tmin=(1272.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_S) + radical(Allyl_P) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH][C]([CH2])[CH]O[O](28518)',
    structure = SMILES('[CH][CH]C=C([CH2])[CH]O[O]'),
    E0 = (836.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376262,0.0700588,-5.57029e-05,1.80666e-08,-9.67949e-13,100794,31.7007], Tmin=(100,'K'), Tmax=(1071.59,'K')), NASAPolynomial(coeffs=[17.2536,0.0217563,-8.6618e-06,1.59954e-09,-1.12117e-13,96332.9,-54.8267], Tmin=(1071.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(836.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_S) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]1[CH]OOC(=C)[CH]1(28519)',
    structure = SMILES('[CH2][C]1[CH]C(=C)[CH]OO1'),
    E0 = (525.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889376,0.0461876,3.31248e-05,-7.669e-08,3.1241e-11,63332.1,24.4733], Tmin=(100,'K'), Tmax=(1028.83,'K')), NASAPolynomial(coeffs=[20.3024,0.023814,-1.16768e-05,2.50856e-09,-1.9428e-13,56527.1,-83.3937], Tmin=(1028.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(C2CsJOOC) + radical(C=CCJCO) + radical(C=CCJO) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][C]1[CH]OOC[C]=C1(28520)',
    structure = SMILES('[CH2]C1=C[C]COO[CH]1'),
    E0 = (633.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17802,0.00257979,0.000226578,-3.42547e-07,1.45245e-10,76395.4,26.9266], Tmin=(100,'K'), Tmax=(908.092,'K')), NASAPolynomial(coeffs=[43.7437,-0.0294941,2.28317e-05,-4.49467e-09,2.90879e-13,62256.5,-209.608], Tmin=(908.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(C=CCJO) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([C]=C=C)[CH]O[O](28521)',
    structure = SMILES('[CH2]C#CC([CH2])[CH]O[O]'),
    E0 = (688.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2100,2250,500,550,180,1768.02],'cm^-1')),
        HinderedRotor(inertia=(0.00177735,'amu*angstrom^2'), symmetry=1, barrier=(3.94281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82791,'amu*angstrom^2'), symmetry=1, barrier=(65.0192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82632,'amu*angstrom^2'), symmetry=1, barrier=(64.9826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171429,'amu*angstrom^2'), symmetry=1, barrier=(3.94148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.82731,'amu*angstrom^2'), symmetry=1, barrier=(65.0055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510862,0.0812233,-0.000112059,8.77399e-08,-2.71288e-11,82898.6,33.0421], Tmin=(100,'K'), Tmax=(906.845,'K')), NASAPolynomial(coeffs=[9.60304,0.0310148,-1.22972e-05,2.11388e-09,-1.36269e-13,81665,-7.64273], Tmin=(906.845,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=CC([CH2])[CH]O[O](28522)',
    structure = SMILES('C#C[CH]C([CH2])[CH]O[O]'),
    E0 = (696.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,1551.27,1551.91],'cm^-1')),
        HinderedRotor(inertia=(0.247719,'amu*angstrom^2'), symmetry=1, barrier=(5.69555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79366,'amu*angstrom^2'), symmetry=1, barrier=(64.2318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248007,'amu*angstrom^2'), symmetry=1, barrier=(5.70217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79528,'amu*angstrom^2'), symmetry=1, barrier=(64.269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79146,'amu*angstrom^2'), symmetry=1, barrier=(64.1813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.144886,0.0901956,-0.000133559,1.05894e-07,-3.21264e-11,83907.8,31.9594], Tmin=(100,'K'), Tmax=(956.695,'K')), NASAPolynomial(coeffs=[10.9734,0.028192,-1.01147e-05,1.59528e-09,-9.52762e-14,82601.5,-15.8017], Tmin=(956.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(Sec_Propargyl) + radical(CCsJOOH) + radical(Isobutyl)"""),
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
    E0 = (589.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (979.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (753.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1139.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1118.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1181.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1046.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (987.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (531.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (531.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (531.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (531.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (586.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (601.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (766.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (762.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (720.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (730.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (744.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (760.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (726.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (545.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (895.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (919.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (934.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (982.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (738.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (739.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (752.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (760.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (717.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (684.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (603.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (770.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (754.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1023.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (719.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (661.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (639.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (660.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (789.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (889.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (798.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (1051.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (804.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (837.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (976.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (725.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (657.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (547.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (914.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (766.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (970.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (530.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (531.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (615.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (578.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (661.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (639.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (696.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (825.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (778.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (1000.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (750.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (564.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (530.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (658.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (896.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (1187.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (854.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (859.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (586.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (666.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (865.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (1011.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['C3H3(5450)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(65.4457,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 65.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH2]C([CH][C]=C)=C[O](25046)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O2(2)', '[CH]=C([CH2])C=[C][CH2](16906)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]O[O](21387)', '[CH2][C]=C[C]=C(17191)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', 'C=[C]C=[C][CH]O[O](28485)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H2CC(T)(1341)', '[CH]C(=C)[CH]O[O](22675)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2]C([C]O[O])=C[C]=C(28486)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C([CH]O[O])=C[C]=C(20221)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]=CC(=C)C1OO1(23665)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['C=[C]C=C1CC1O[O](28487)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1=CC(=C)C1O[O](22695)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_H/NonDeO;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1=CC(=CO[O])C1(22672)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C(=C=C=C)CO[O](28488)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['C=C=C=C(C)[CH]O[O](28489)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]1C([C]=C)C1O[O](28490)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(242.814,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 239.7 to 242.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['C=[C]C1C[C]1[CH]O[O](28491)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(238.957,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 236.2 to 239.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1=CC1([CH2])[CH]O[O](22938)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(197.301,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 197.2 to 197.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]1[CH]OOC1[C]=C(28492)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(206.715,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 202.7 to 206.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1([CH][C]=C)[CH]OO1(28493)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(221.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]C(=C=C=C)[CH]O[O](28494)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH]=C=CC(=C)[CH]O[O](22639)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.156e+09,'cm^3/(mol*s)'), n=1.502, Ea=(9.92026,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 194 used for Ct-H_Ct-Cd;HJ
Exact match found for rate rule [Ct-H_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH2]C(C=O)=C[C]=C(28495)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(20.92,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(T)(20)', 'C=[C]C=C=CO[O](28496)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(10.4588,'m^3/(mol*s)'), n=1.99333, Ea=(11.9069,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;YJ] for rate rule [Ca_Cds;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]O[O](21387)', 'C=[C]C=C=C(17211)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.40891,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]C(=[C][C]=C)[CH]O[O](28497)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH]=[C][CH]C(=C)[CH]O[O](22641)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(=[C]C=C)[CH]O[O](28498)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C[CH]C(=C)[CH]O[O](22642)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=[C]C(=C)CO[O](20234)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C][C]=C(C)[CH]O[O](28499)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C][CH]C(=C)CO[O](20241)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.1728e+06,'s^-1'), n=1.70245, Ea=(63.8935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cs_H_out_1H] + [R5Hall;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=[C][CH]C(C)=CO[O](22653)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=[C][C]=C)[CH]OO(28500)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C][CH]C(=C)[CH]OO(22643)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C]=[C]C([CH2])[CH]O[O](28501)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]=C[C]1CC1O[O](28502)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_csHO]
Euclidian distance = 2.8284271247461903
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1=C[C]([CH]O[O])C1(28503)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[O]O[CH][C]1C=[C]CC1(28504)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.49459e+10,'s^-1'), n=0.314867, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]=C[C]1[CH]OOC1(28505)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.20877e+14,'s^-1'), n=-0.684234, Ea=(136.856,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secDe_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1([CH]O[O])C=[C]C1(22637)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.36716e+08,'s^-1'), n=1.01412, Ea=(266.197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 265.8 to 266.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C][CH2](16918)', 'C=C=CO[O](16806)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(0.00267131,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C3H3(5450)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C][CH2](16918)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C=[C][CH2])CO[O](19510)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C([C]=[C]C)[CH]O[O](28506)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]C([CH]O[O])=CC=C(20227)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C([CH]O[O])C=[C]C(20236)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]C([CH]OO)=C[C]=C(20226)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['C=[C]C1C(=C)C1O[O](22556)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(23.7316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C([CH]O[O])C[C]=C(20080)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C=CC[C][CH]O[O](20294)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]=CC1=COOC1(19223)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1=COOC([CH2])=C1(28507)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['O(4)', '[CH2][C]=CC1=COC1(20232)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(92.0547,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['O(4)', '[CH2]C1=COC([CH2])=C1(28508)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_DSD;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]1C=C([CH2])C1O[O](28509)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]1C=[C]CC1O[O](28510)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(4.49459e+10,'s^-1'), n=0.314867, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]=C[C]([CH2])C1OO1(28511)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(173.289,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 171.5 to 173.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2][C]=CC([CH2])=[C]OO(28512)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]=CC(C)=[C]O[O](28513)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C=CC([CH2])=[C]O[O](28514)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C(=[C]O[O])C=[C]C(28515)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['C=[C]C1CC1=CO[O](28516)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(40.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2]C1=COOC1[C]=C(22696)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['O(4)', '[CH2]C1=COC1[C]=C(22697)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(134.889,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO] for rate rule [R3OO_SD;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction68',
    reactants = ['C=[C]CC(=C)[C]O[O](22313)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['C3H3(5450)', '[CH2][C][CH]O[O](22548)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]C=[C][C]([CH2])[CH]O[O](28517)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]=C[CH][C]([CH2])[CH]O[O](28518)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]1[CH]OOC(=C)[CH]1(28519)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    products = ['[CH2][C]1[CH]OOC[C]=C1(28520)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(3.22e+12,'s^-1'), n=-0.622, Ea=(142.884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_CdCdd;radadd_intra] for rate rule [R7;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH2]C([C]=C=C)[CH]O[O](28521)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]=C=CC([CH2])[CH]O[O](28522)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #5017',
    isomers = [
        '[CH2][C]=CC(=C)[CH]O[O](20233)',
    ],
    reactants = [
        ('C3H3(5450)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #5017',
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


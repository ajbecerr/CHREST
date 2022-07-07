species(
    label = '[CH2]C(=C)C(=C)[CH]O[O](20115)',
    structure = SMILES('[CH2]C(=C)C(=C)[CH]O[O]'),
    E0 = (357.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,410.701,411.05],'cm^-1')),
        HinderedRotor(inertia=(0.197226,'amu*angstrom^2'), symmetry=1, barrier=(23.7177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197435,'amu*angstrom^2'), symmetry=1, barrier=(23.7194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198414,'amu*angstrom^2'), symmetry=1, barrier=(23.7118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282394,'amu*angstrom^2'), symmetry=1, barrier=(33.9178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.144508,0.0724169,-4.57391e-05,-1.62991e-09,8.31116e-12,43196.4,28.4949], Tmin=(100,'K'), Tmax=(970.222,'K')), NASAPolynomial(coeffs=[19.9217,0.0189297,-6.4114e-06,1.14464e-09,-8.18386e-14,38038.5,-73.1243], Tmin=(970.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=C)C([C]=C)O[O](21270)',
    structure = SMILES('[CH2]C(=C)C([C]=C)O[O]'),
    E0 = (484.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,970.938],'cm^-1')),
        HinderedRotor(inertia=(0.550969,'amu*angstrom^2'), symmetry=1, barrier=(12.6679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0405387,'amu*angstrom^2'), symmetry=1, barrier=(12.6707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551035,'amu*angstrom^2'), symmetry=1, barrier=(12.6694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29794,'amu*angstrom^2'), symmetry=1, barrier=(29.8421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3959.8,'J/mol'), sigma=(6.603,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=618.51 K, Pc=31.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661558,0.0755181,-7.52483e-05,4.0531e-08,-8.99755e-12,58415.9,30.2303], Tmin=(100,'K'), Tmax=(1071.96,'K')), NASAPolynomial(coeffs=[12.2832,0.0321531,-1.45686e-05,2.79431e-09,-1.96862e-13,55924.3,-26.6452], Tmin=(1071.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(=C)[CH]O[O](20082)',
    structure = SMILES('C=[C]CC(=C)[CH]O[O]'),
    E0 = (475.206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,771.66,771.665],'cm^-1')),
        HinderedRotor(inertia=(0.179678,'amu*angstrom^2'), symmetry=1, barrier=(4.13115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0362529,'amu*angstrom^2'), symmetry=1, barrier=(15.3189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.666273,'amu*angstrom^2'), symmetry=1, barrier=(15.3189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39634,'amu*angstrom^2'), symmetry=1, barrier=(32.1047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3938.88,'J/mol'), sigma=(6.5717,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.24 K, Pc=31.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544713,0.0701626,-5.99611e-05,2.64267e-08,-4.69146e-12,57283.2,32.3323], Tmin=(100,'K'), Tmax=(1342.53,'K')), NASAPolynomial(coeffs=[15.3522,0.0260448,-1.0669e-05,1.94969e-09,-1.33493e-13,53307.3,-43.4669], Tmin=(1342.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
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
    label = '[CH2]C(=C)C(=C)[CH][O](23578)',
    structure = SMILES('[CH2]C(=C)C([CH2])=C[O]'),
    E0 = (254.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43333,'amu*angstrom^2'), symmetry=1, barrier=(32.9551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4326,'amu*angstrom^2'), symmetry=1, barrier=(32.9383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43095,'amu*angstrom^2'), symmetry=1, barrier=(32.9003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40155,0.0603326,-6.20328e-06,-5.23668e-08,2.94114e-11,30745.2,23.9436], Tmin=(100,'K'), Tmax=(925.932,'K')), NASAPolynomial(coeffs=[23.7407,0.00806211,-1.82938e-07,-6.85501e-11,2.11033e-16,24341.7,-98.098], Tmin=(925.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH]C(=C)C([CH2])=C(16082)',
    structure = SMILES('[CH]C(=C)C([CH2])=C'),
    E0 = (540.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,195.967,195.967,195.967,195.967],'cm^-1')),
        HinderedRotor(inertia=(1.80853,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80854,'amu*angstrom^2'), symmetry=1, barrier=(49.2857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3249.64,'J/mol'), sigma=(5.69809,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.59 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01927,0.0540672,-1.36981e-05,-2.39878e-08,1.40944e-11,65175.7,20.5913], Tmin=(100,'K'), Tmax=(963.047,'K')), NASAPolynomial(coeffs=[14.94,0.0233812,-8.1649e-06,1.43761e-09,-1.00383e-13,61236.2,-52.5768], Tmin=(963.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C(=C)[C]=C(15786)',
    structure = SMILES('[CH2]C(=C)[C]=C'),
    E0 = (406.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.68652,'amu*angstrom^2'), symmetry=1, barrier=(38.7765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67506,'amu*angstrom^2'), symmetry=1, barrier=(38.5129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7794,0.0407457,-1.28305e-05,-1.73438e-08,1.13623e-11,49017.2,16.4439], Tmin=(100,'K'), Tmax=(926.058,'K')), NASAPolynomial(coeffs=[12.6514,0.0145843,-4.14448e-06,6.55926e-10,-4.4289e-14,46111.7,-39.9871], Tmin=(926.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
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
    label = '[CH2]C(=C)C([CH2])=[C]O[O](22952)',
    structure = SMILES('[CH2]C(=C)C([CH2])=[C]O[O]'),
    E0 = (629.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.619135,'amu*angstrom^2'), symmetry=1, barrier=(14.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.619074,'amu*angstrom^2'), symmetry=1, barrier=(14.2337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0547146,'amu*angstrom^2'), symmetry=1, barrier=(14.2397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310585,'amu*angstrom^2'), symmetry=1, barrier=(80.6864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.178991,0.080722,-9.01424e-05,5.17233e-08,-1.16538e-11,75814.5,30.5694], Tmin=(100,'K'), Tmax=(1088.31,'K')), NASAPolynomial(coeffs=[16.3954,0.0211199,-7.99371e-06,1.40138e-09,-9.41305e-14,72284.8,-49.0375], Tmin=(1088.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=C)C(=C)[CH]O[O](20338)',
    structure = SMILES('[CH]C(=C)C(=C)[CH]O[O]'),
    E0 = (577.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3958.16,'J/mol'), sigma=(6.59471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=618.26 K, Pc=31.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123104,0.0750846,-5.47386e-05,1.21598e-08,2.17045e-12,69556.9,30.0934], Tmin=(100,'K'), Tmax=(1003.87,'K')), NASAPolynomial(coeffs=[17.7075,0.0248468,-9.30081e-06,1.66062e-09,-1.14836e-13,65027.2,-59.7856], Tmin=(1003.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)C(=C)C1OO1(23579)',
    structure = SMILES('[CH2]C(=C)C(=C)C1OO1'),
    E0 = (204.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333113,0.0628675,-8.34469e-06,-5.00916e-08,2.86951e-11,24731.2,25.358], Tmin=(100,'K'), Tmax=(918.353,'K')), NASAPolynomial(coeffs=[22.7151,0.0119829,-1.35061e-06,8.86165e-11,-7.85401e-15,18655.1,-91.415], Tmin=(918.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(dioxirane) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1CC(O[O])C1=C(22292)',
    structure = SMILES('C=C1CC(O[O])C1=C'),
    E0 = (236.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901638,0.0548963,-1.04432e-05,-2.74632e-08,1.4705e-11,28611.1,24.4924], Tmin=(100,'K'), Tmax=(998.324,'K')), NASAPolynomial(coeffs=[16.0547,0.0241158,-9.17026e-06,1.72071e-09,-1.24266e-13,24093.9,-56.0576], Tmin=(998.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C([CH]O[O])C([CH2])[CH2](20017)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])[CH2]'),
    E0 = (761.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,479.555],'cm^-1')),
        HinderedRotor(inertia=(0.0976643,'amu*angstrom^2'), symmetry=1, barrier=(15.9403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104629,'amu*angstrom^2'), symmetry=1, barrier=(80.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104631,'amu*angstrom^2'), symmetry=1, barrier=(80.6188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494045,'amu*angstrom^2'), symmetry=1, barrier=(80.6171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.50634,'amu*angstrom^2'), symmetry=1, barrier=(80.6177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269947,0.0772992,-7.47279e-05,3.4959e-08,-5.25048e-12,91760.2,35.3693], Tmin=(100,'K'), Tmax=(923.393,'K')), NASAPolynomial(coeffs=[15.5168,0.023695,-7.86367e-06,1.27778e-09,-8.25316e-14,88414,-39.8451], Tmin=(923.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(761.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C([CH2])[CH]O[O](18233)',
    structure = SMILES('[CH]C(=C)C([CH2])[CH]O[O]'),
    E0 = (744.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156733,0.0884539,-0.000107567,7.84661e-08,-2.37232e-11,89723.9,33.8976], Tmin=(100,'K'), Tmax=(821.235,'K')), NASAPolynomial(coeffs=[9.82831,0.040098,-1.69643e-05,3.06474e-09,-2.06049e-13,88177.4,-10.6011], Tmin=(821.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH]O[O])[C]([CH2])C(20107)',
    structure = SMILES('[CH]C([CH]O[O])=C([CH2])C'),
    E0 = (601.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.22062,0.0780973,-6.51285e-05,2.85855e-08,-5.13184e-12,72524.8,31.3143], Tmin=(100,'K'), Tmax=(1315.28,'K')), NASAPolynomial(coeffs=[15.0757,0.0329204,-1.36068e-05,2.4711e-09,-1.68169e-13,68617.1,-44.4238], Tmin=(1315.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])[C](C)[CH]O[O](23580)',
    structure = SMILES('[CH]C([CH2])=C(C)[CH]O[O]'),
    E0 = (601.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.22062,0.0780973,-6.51285e-05,2.85855e-08,-5.13184e-12,72524.8,31.3143], Tmin=(100,'K'), Tmax=(1315.28,'K')), NASAPolynomial(coeffs=[15.0757,0.0329204,-1.36068e-05,2.4711e-09,-1.68169e-13,68617.1,-44.4238], Tmin=(1315.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)[C]1CC1O[O](23581)',
    structure = SMILES('[CH2]C(=C)[C]1CC1O[O]'),
    E0 = (430.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874776,0.053502,-4.42175e-07,-4.1827e-08,2.08347e-11,51888.9,27.9594], Tmin=(100,'K'), Tmax=(973.075,'K')), NASAPolynomial(coeffs=[17.3305,0.0220625,-7.78807e-06,1.44202e-09,-1.05415e-13,46972.3,-59.7884], Tmin=(973.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([CH]O[O])[C]1CC1(23582)',
    structure = SMILES('C=C([CH]O[O])[C]1CC1'),
    E0 = (394.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19652,0.0440152,2.5944e-05,-6.90359e-08,3.06128e-11,47614.5,26.5762], Tmin=(100,'K'), Tmax=(958.192,'K')), NASAPolynomial(coeffs=[16.9906,0.0215808,-7.031e-06,1.28381e-09,-9.52377e-14,42590.9,-59.3667], Tmin=(958.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(Allyl_T) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C]1CC(O[O])C1=C(23583)',
    structure = SMILES('[CH2]C1=C([CH2])C(O[O])C1'),
    E0 = (420.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.592334,0.0595603,-1.19585e-05,-3.3296e-08,1.84788e-11,50657.9,27.7008], Tmin=(100,'K'), Tmax=(975.857,'K')), NASAPolynomial(coeffs=[18.9358,0.0204337,-7.24902e-06,1.35587e-09,-9.99948e-14,45360.6,-69.1457], Tmin=(975.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cyclobutene) + radical(ROOJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1CC[C]1[CH]O[O](23584)',
    structure = SMILES('[CH2]C1=C([CH]O[O])CC1'),
    E0 = (405.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864613,0.0519717,8.13531e-06,-5.26809e-08,2.50106e-11,48901.6,28.2211], Tmin=(100,'K'), Tmax=(970.332,'K')), NASAPolynomial(coeffs=[18.2605,0.0211325,-7.37354e-06,1.38376e-09,-1.02975e-13,43601.6,-65.0948], Tmin=(970.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cyclobutene) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)[C]1[CH]OOC1(23585)',
    structure = SMILES('[CH2]C([CH2])=C1[CH]OOC1'),
    E0 = (323.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33487,0.0331543,7.26247e-05,-1.228e-07,5.02268e-11,39049.9,26.8281], Tmin=(100,'K'), Tmax=(959.822,'K')), NASAPolynomial(coeffs=[19.7741,0.0200233,-6.42479e-06,1.26429e-09,-1.01154e-13,32575.4,-76.6625], Tmin=(959.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(Cyclopentane) + radical(C=CCJO) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]1COO[CH]C1=C(23586)',
    structure = SMILES('[CH2]C1=C([CH2])COO[CH]1'),
    E0 = (310.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10624,0.0373315,6.57346e-05,-1.20988e-07,5.1052e-11,37489.8,22.8566], Tmin=(100,'K'), Tmax=(952.04,'K')), NASAPolynomial(coeffs=[21.9943,0.0156241,-4.13623e-06,8.16072e-10,-6.99218e-14,30519,-92.6113], Tmin=(952.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(36dihydro12dioxin) + radical(C=CCJO) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([CH2])C(=C)C1O[O](22358)',
    structure = SMILES('[CH2]C1([CH2])C(=C)C1O[O]'),
    E0 = (583.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.34476,0.0699859,-5.15852e-05,1.30417e-08,8.25606e-13,70281.6,29.6146], Tmin=(100,'K'), Tmax=(1062.79,'K')), NASAPolynomial(coeffs=[17.61,0.0223835,-8.92735e-06,1.66873e-09,-1.18192e-13,65630.2,-59.349], Tmin=(1062.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C1([CH]O[O])CC1=C(22319)',
    structure = SMILES('[CH2]C1([CH]O[O])CC1=C'),
    E0 = (579.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510249,0.0727039,-6.4046e-05,2.90978e-08,-5.33983e-12,69856.6,27.5842], Tmin=(100,'K'), Tmax=(1297.09,'K')), NASAPolynomial(coeffs=[15.1982,0.0274086,-1.16647e-05,2.17514e-09,-1.50754e-13,66046.3,-47.0973], Tmin=(1297.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(CCsJOOH) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(=C)C1([CH2])[CH]OO1(23587)',
    structure = SMILES('[CH2]C(=C)C1([CH2])[CH]OO1'),
    E0 = (516.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275456,0.0755164,-6.89684e-05,3.27719e-08,-6.22808e-12,62267.9,27.0675], Tmin=(100,'K'), Tmax=(1268.39,'K')), NASAPolynomial(coeffs=[16.0819,0.0256688,-1.00185e-05,1.78771e-09,-1.21064e-13,58258.1,-52.9475], Tmin=(1268.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(CJCOOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([CH2])OO[CH]C1=C(23588)',
    structure = SMILES('[CH2]C1([CH2])OO[CH]C1=C'),
    E0 = (432.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02549,0.0384364,6.54799e-05,-1.24588e-07,5.36318e-11,52156.9,31.0643], Tmin=(100,'K'), Tmax=(941.371,'K')), NASAPolynomial(coeffs=[23.2187,0.0124828,-2.07224e-06,3.77735e-10,-3.8088e-14,44950,-90.7493], Tmin=(941.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJO) + radical(CJCOOH) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C(=C)C(=C)C=O(23589)',
    structure = SMILES('[CH2]C(=C)C(=C)C=O'),
    E0 = (93.8763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.09006,'amu*angstrom^2'), symmetry=1, barrier=(25.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09235,'amu*angstrom^2'), symmetry=1, barrier=(25.1153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09357,'amu*angstrom^2'), symmetry=1, barrier=(25.1434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544084,0.0637848,-3.72225e-05,-3.29144e-09,7.01794e-12,11425.8,21.827], Tmin=(100,'K'), Tmax=(1018.03,'K')), NASAPolynomial(coeffs=[18.5295,0.0182982,-7.30345e-06,1.41242e-09,-1.03727e-13,6459.03,-71.6726], Tmin=(1018.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(93.8763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P)"""),
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
    label = '[CH]=C([CH]O[O])C([CH2])=C(20113)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])=C'),
    E0 = (605.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,349.682],'cm^-1')),
        HinderedRotor(inertia=(0.282932,'amu*angstrom^2'), symmetry=1, barrier=(24.3121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393143,'amu*angstrom^2'), symmetry=1, barrier=(33.8522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281759,'amu*angstrom^2'), symmetry=1, barrier=(24.3754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282006,'amu*angstrom^2'), symmetry=1, barrier=(24.3557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0702979,0.0759462,-6.24266e-05,1.58034e-08,2.36754e-12,72916.1,29.1644], Tmin=(100,'K'), Tmax=(970.967,'K')), NASAPolynomial(coeffs=[20.1627,0.0160951,-5.37621e-06,9.45567e-10,-6.69437e-14,67933.8,-72.7421], Tmin=(970.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(CO[O])C([CH2])=C(19489)',
    structure = SMILES('[CH]=C(CO[O])C([CH2])=C'),
    E0 = (487.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,222.496],'cm^-1')),
        HinderedRotor(inertia=(0.424134,'amu*angstrom^2'), symmetry=1, barrier=(14.8406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423787,'amu*angstrom^2'), symmetry=1, barrier=(14.8301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419331,'amu*angstrom^2'), symmetry=1, barrier=(14.8292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.37592,'amu*angstrom^2'), symmetry=1, barrier=(116.708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3958.16,'J/mol'), sigma=(6.59471,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=618.26 K, Pc=31.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102827,0.0841266,-9.21653e-05,5.2826e-08,-1.20325e-11,58799.6,28.4536], Tmin=(100,'K'), Tmax=(1069.76,'K')), NASAPolynomial(coeffs=[15.6908,0.0258405,-1.04367e-05,1.89293e-09,-1.29462e-13,55464.6,-47.8002], Tmin=(1069.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C)C(=C)[CH]O[O](23590)',
    structure = SMILES('[CH]=C(C)C(=C)[CH]O[O]'),
    E0 = (453.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,272.974],'cm^-1')),
        HinderedRotor(inertia=(0.337125,'amu*angstrom^2'), symmetry=1, barrier=(17.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337125,'amu*angstrom^2'), symmetry=1, barrier=(17.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337125,'amu*angstrom^2'), symmetry=1, barrier=(17.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.520285,'amu*angstrom^2'), symmetry=1, barrier=(27.5114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.135861,0.0807965,-7.88232e-05,3.87154e-08,-7.43059e-12,54701.6,30.3215], Tmin=(100,'K'), Tmax=(1275.77,'K')), NASAPolynomial(coeffs=[19.335,0.0197477,-7.04398e-06,1.20624e-09,-8.02339e-14,49733.5,-68.3564], Tmin=(1275.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C(=C)C(20117)',
    structure = SMILES('[CH]=C([CH]O[O])C(=C)C'),
    E0 = (453.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,272.974],'cm^-1')),
        HinderedRotor(inertia=(0.337125,'amu*angstrom^2'), symmetry=1, barrier=(17.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337125,'amu*angstrom^2'), symmetry=1, barrier=(17.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.337125,'amu*angstrom^2'), symmetry=1, barrier=(17.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.520285,'amu*angstrom^2'), symmetry=1, barrier=(27.5114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.135861,0.0807965,-7.88232e-05,3.87154e-08,-7.43059e-12,54701.6,30.3215], Tmin=(100,'K'), Tmax=(1275.77,'K')), NASAPolynomial(coeffs=[19.335,0.0197477,-7.04398e-06,1.20624e-09,-8.02339e-14,49733.5,-68.3564], Tmin=(1275.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C(=C)CO[O](18231)',
    structure = SMILES('[CH]C(=C)C(=C)CO[O]'),
    E0 = (459.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.27257,0.0818734,-7.95608e-05,4.27495e-08,-9.46494e-12,55435.3,28.964], Tmin=(100,'K'), Tmax=(1079,'K')), NASAPolynomial(coeffs=[12.8514,0.0352416,-1.47341e-05,2.69574e-09,-1.84613e-13,52720.8,-32.6781], Tmin=(1079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH]OO)C([CH2])=C(20114)',
    structure = SMILES('[CH]=C([CH]OO)C([CH2])=C'),
    E0 = (453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,320.726],'cm^-1')),
        HinderedRotor(inertia=(0.422028,'amu*angstrom^2'), symmetry=1, barrier=(30.8,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421971,'amu*angstrom^2'), symmetry=1, barrier=(30.7997,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421982,'amu*angstrom^2'), symmetry=1, barrier=(30.8002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421951,'amu*angstrom^2'), symmetry=1, barrier=(30.8001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421976,'amu*angstrom^2'), symmetry=1, barrier=(30.8001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.220596,0.0798449,-5.87298e-05,7.68679e-09,5.67701e-12,54646.8,29.45], Tmin=(100,'K'), Tmax=(981.833,'K')), NASAPolynomial(coeffs=[22.052,0.0172408,-6.06909e-06,1.11561e-09,-8.11928e-14,48917.1,-84.5], Tmin=(981.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH2])C(=C)[CH]OO(23591)',
    structure = SMILES('[CH]C(=C)C(=C)[CH]OO'),
    E0 = (425.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162548,0.0789167,-5.07853e-05,3.6744e-09,5.65416e-12,51287.3,30.3605], Tmin=(100,'K'), Tmax=(1007.44,'K')), NASAPolynomial(coeffs=[19.5932,0.0260002,-9.99868e-06,1.83195e-09,-1.29198e-13,46011.6,-71.5243], Tmin=(1007.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)C[C]=CO[O](21271)',
    structure = SMILES('[CH2]C(=C)C[C]=CO[O]'),
    E0 = (506.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,232.746,233.291],'cm^-1')),
        HinderedRotor(inertia=(0.375271,'amu*angstrom^2'), symmetry=1, barrier=(14.5869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00306672,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378705,'amu*angstrom^2'), symmetry=1, barrier=(14.5889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374977,'amu*angstrom^2'), symmetry=1, barrier=(14.5877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4016.09,'J/mol'), sigma=(6.63908,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=627.30 K, Pc=31.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40762,0.0763236,-7.65235e-05,4.04687e-08,-8.57695e-12,61078.7,31.1977], Tmin=(100,'K'), Tmax=(1141.31,'K')), NASAPolynomial(coeffs=[14.7396,0.0260945,-1.05096e-05,1.90914e-09,-1.3076e-13,57807.2,-39.8406], Tmin=(1141.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(=C)[C]=CO[O](23592)',
    structure = SMILES('[CH2]C(=C)[C]=CO[O]'),
    E0 = (474.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.10879,'amu*angstrom^2'), symmetry=1, barrier=(25.4933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10842,'amu*angstrom^2'), symmetry=1, barrier=(25.4848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10884,'amu*angstrom^2'), symmetry=1, barrier=(25.4945,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696334,0.0684485,-7.19061e-05,3.41268e-08,-4.68206e-12,57196,24.9934], Tmin=(100,'K'), Tmax=(896.38,'K')), NASAPolynomial(coeffs=[15.8292,0.0144577,-4.21294e-06,6.30413e-10,-3.91561e-14,53939.1,-49.3926], Tmin=(896.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1CCC1=CO[O](22280)',
    structure = SMILES('C=C1CCC1=CO[O]'),
    E0 = (253.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932591,0.0546674,-1.08985e-05,-2.80308e-08,1.55159e-11,30654.8,25.6395], Tmin=(100,'K'), Tmax=(975.027,'K')), NASAPolynomial(coeffs=[15.9472,0.0229489,-8.06735e-06,1.46188e-09,-1.04469e-13,26306.7,-53.7009], Tmin=(975.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(=C)C1=COOC1(19219)',
    structure = SMILES('[CH2]C(=C)C1=COOC1'),
    E0 = (169.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.66913,0.0545619,8.09924e-06,-5.90762e-08,2.90518e-11,20555.1,25.0771], Tmin=(100,'K'), Tmax=(951.985,'K')), NASAPolynomial(coeffs=[20.7067,0.0165392,-4.73888e-06,8.59531e-10,-6.65942e-14,14647.9,-81.5952], Tmin=(951.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolene) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1=COOCC1=C(23593)',
    structure = SMILES('[CH2]C1=COOCC1=C'),
    E0 = (152.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07653,0.039952,5.69247e-05,-1.04696e-07,4.25926e-11,18476.6,19.3583], Tmin=(100,'K'), Tmax=(984.567,'K')), NASAPolynomial(coeffs=[19.812,0.0236581,-9.39182e-06,1.92096e-09,-1.49404e-13,11887.8,-85.4628], Tmin=(984.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1=COC1(20105)',
    structure = SMILES('[CH2]C(=C)C1=COC1'),
    E0 = (168.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817047,0.0428577,5.18198e-05,-1.18647e-07,5.4542e-11,20458,19.089], Tmin=(100,'K'), Tmax=(921.302,'K')), NASAPolynomial(coeffs=[26.4264,0.00208398,3.56159e-06,-7.69106e-10,4.44596e-14,12750.8,-118.58], Tmin=(921.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1=COCC1=C(23594)',
    structure = SMILES('[CH2]C1=COCC1=C'),
    E0 = (73.2461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46733,0.0213591,0.000115654,-1.86847e-07,7.95507e-11,8932.64,18.2167], Tmin=(100,'K'), Tmax=(917.053,'K')), NASAPolynomial(coeffs=[26.9827,-0.000818946,6.1667e-06,-1.2896e-09,7.80805e-14,505.631,-123.102], Tmin=(917.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.2461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)[C]O[O](22189)',
    structure = SMILES('[CH2]C(=[C]O[O])C([CH2])[CH2]'),
    E0 = (785.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,1457.81],'cm^-1')),
        HinderedRotor(inertia=(0.199995,'amu*angstrom^2'), symmetry=1, barrier=(5.59937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212164,'amu*angstrom^2'), symmetry=1, barrier=(5.60775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20491,'amu*angstrom^2'), symmetry=1, barrier=(5.60623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00374311,'amu*angstrom^2'), symmetry=1, barrier=(5.65373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41189,'amu*angstrom^2'), symmetry=1, barrier=(64.758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625155,0.078529,-8.63558e-05,4.29342e-08,-3.05144e-12,94648.1,35.9269], Tmin=(100,'K'), Tmax=(679.072,'K')), NASAPolynomial(coeffs=[11.6406,0.0289451,-1.06287e-05,1.77161e-09,-1.13075e-13,92799.3,-15.5505], Tmin=(679.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C)C([CH2])=[C]O[O](23595)',
    structure = SMILES('[CH2][C](C)C([CH2])=[C]O[O]'),
    E0 = (712.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2099.66],'cm^-1')),
        HinderedRotor(inertia=(0.0218721,'amu*angstrom^2'), symmetry=1, barrier=(68.4273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245264,'amu*angstrom^2'), symmetry=1, barrier=(6.08391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00482207,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24527,'amu*angstrom^2'), symmetry=1, barrier=(6.08411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.75891,'amu*angstrom^2'), symmetry=1, barrier=(68.4273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654896,0.0747715,-8.36093e-05,5.44299e-08,-1.45079e-11,85857.6,33.6429], Tmin=(100,'K'), Tmax=(910.104,'K')), NASAPolynomial(coeffs=[10.511,0.0314522,-1.22108e-05,2.12828e-09,-1.40685e-13,84063.6,-12.9784], Tmin=(910.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(712.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Allyl_T) + radical(Isobutyl) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=C)[C]([CH2])CO[O](18232)',
    structure = SMILES('[CH]C([CH2])=C([CH2])CO[O]'),
    E0 = (636.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.3724,0.0823155,-8.10295e-05,4.5409e-08,-1.06925e-11,76627,30.4615], Tmin=(100,'K'), Tmax=(1007.15,'K')), NASAPolynomial(coeffs=[11.144,0.0395355,-1.73157e-05,3.23516e-09,-2.23987e-13,74457.3,-21.5822], Tmin=(1007.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C](C([CH2])=C)C1OO1(23596)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])C1OO1'),
    E0 = (380.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.306156,0.0648962,-1.57724e-05,-3.91506e-08,2.3701e-11,45928.3,27.9971], Tmin=(100,'K'), Tmax=(928.703,'K')), NASAPolynomial(coeffs=[21.5179,0.0153896,-3.41234e-06,5.03795e-10,-3.68312e-14,40183.5,-82.4858], Tmin=(928.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + ring(dioxirane) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([CH2])CC1=CO[O](22341)',
    structure = SMILES('[CH2]C1([CH2])CC1=CO[O]'),
    E0 = (600.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390457,0.0696072,-5.16418e-05,1.21619e-08,1.67067e-12,72324.6,30.7073], Tmin=(100,'K'), Tmax=(1016.26,'K')), NASAPolynomial(coeffs=[17.3031,0.0215411,-8.00566e-06,1.4517e-09,-1.01798e-13,67931.7,-55.8599], Tmin=(1016.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])=[C]OO(23597)',
    structure = SMILES('[CH2]C(=C)C([CH2])=[C]OO'),
    E0 = (477.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.48689,'amu*angstrom^2'), symmetry=1, barrier=(34.1866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0679271,'amu*angstrom^2'), symmetry=1, barrier=(18.5522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205275,'amu*angstrom^2'), symmetry=1, barrier=(18.553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.67443,'amu*angstrom^2'), symmetry=1, barrier=(84.4825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0678885,'amu*angstrom^2'), symmetry=1, barrier=(18.5518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.179605,0.0854118,-8.91474e-05,4.69632e-08,-9.69259e-12,57548.1,31.0981], Tmin=(100,'K'), Tmax=(1185.77,'K')), NASAPolynomial(coeffs=[18.7865,0.0214331,-8.21481e-06,1.4614e-09,-9.93399e-14,53050.2,-63.6342], Tmin=(1185.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=C)C(C)=[C]O[O](23598)',
    structure = SMILES('[CH2]C(=C)C(C)=[C]O[O]'),
    E0 = (477.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333506,0.0813346,-9.18665e-05,5.59594e-08,-1.36862e-11,57584.4,29.7394], Tmin=(100,'K'), Tmax=(994.216,'K')), NASAPolynomial(coeffs=[13.5654,0.028099,-1.15484e-05,2.10248e-09,-1.43625e-13,54953.3,-34.0202], Tmin=(994.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C(=C)C(23599)',
    structure = SMILES('[CH2]C(=[C]O[O])C(=C)C'),
    E0 = (477.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,270.457],'cm^-1')),
        HinderedRotor(inertia=(0.234555,'amu*angstrom^2'), symmetry=1, barrier=(12.173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234544,'amu*angstrom^2'), symmetry=1, barrier=(12.1728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23451,'amu*angstrom^2'), symmetry=1, barrier=(12.1729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234528,'amu*angstrom^2'), symmetry=1, barrier=(12.1727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333506,0.0813346,-9.18665e-05,5.59594e-08,-1.36862e-11,57584.4,30.4326], Tmin=(100,'K'), Tmax=(994.216,'K')), NASAPolynomial(coeffs=[13.5654,0.028099,-1.15484e-05,2.10248e-09,-1.43625e-13,54953.3,-33.327], Tmin=(994.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([CH2])C(C)=CO[O](23600)',
    structure = SMILES('[CH]C(=C)C(C)=CO[O]'),
    E0 = (457.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14207,'amu*angstrom^2'), symmetry=1, barrier=(49.2504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14206,'amu*angstrom^2'), symmetry=1, barrier=(49.2502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1421,'amu*angstrom^2'), symmetry=1, barrier=(49.251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14217,'amu*angstrom^2'), symmetry=1, barrier=(49.2527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0999091,0.0847125,-8.29712e-05,4.30593e-08,-8.93439e-12,55132.5,29.6816], Tmin=(100,'K'), Tmax=(1168.9,'K')), NASAPolynomial(coeffs=[16.347,0.0284306,-1.07468e-05,1.86693e-09,-1.24298e-13,51287.5,-52.232], Tmin=(1168.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
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
    E0 = (426.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (579.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (645.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (773.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (587.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (933.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (880.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (840.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (788.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (365.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (366.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (825.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (808.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (626.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (626.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (553.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (451.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (483.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (488.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (494.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (421.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (583.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (579.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (516.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (452.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (357.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (670.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (645.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (832.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (816.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (632.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (582.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (497.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (504.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (486.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (576.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (676.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (912.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (366.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (365.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (366.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (440.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (419.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (849.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (721.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (661.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (404.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (600.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (619.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (613.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (521.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (501.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['C=C=CO[O](16806)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(68.4604,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 68.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=C)C([C]=C)O[O](21270)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]CC(=C)[CH]O[O](20082)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH2]C(=C)C(=C)[CH][O](23578)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O2(2)', '[CH]C(=C)C([CH2])=C(16082)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]O[O](21387)', '[CH2]C(=C)[C]=C(15786)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', 'C=[C]C(=C)[CH]O[O](22489)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2]C(=C)C([CH2])=[C]O[O](22952)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C(=C)C(=C)[CH]O[O](20338)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C(=C)C(=C)C1OO1(23579)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['C=C1CC(O[O])C1=C(22292)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C([CH]O[O])C([CH2])[CH2](20017)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(=C)C([CH2])[CH]O[O](18233)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C([CH]O[O])[C]([CH2])C(20107)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C([CH2])[C](C)[CH]O[O](23580)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C(=C)[C]1CC1O[O](23581)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(195.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_csHO]
Euclidian distance = 2.8284271247461903
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['C=C([CH]O[O])[C]1CC1(23582)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(93.5715,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2][C]1CC(O[O])C1=C(23583)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.3672e+08,'s^-1'), n=0.948854, Ea=(125.32,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['C=C1CC[C]1[CH]O[O](23584)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.43734e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C(=C)[C]1[CH]OOC1(23585)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.20877e+14,'s^-1'), n=-0.684234, Ea=(136.856,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secDe_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2][C]1COO[CH]C1=C(23586)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.9903e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C1([CH2])C(=C)C1O[O](22358)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(225.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 224.6 to 225.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C1([CH]O[O])CC1=C(22319)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.48119e+09,'s^-1'), n=0.736667, Ea=(221.841,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 221.8 to 221.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C(=C)C1([CH2])[CH]OO1(23587)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(158.842,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra_2H;radadd_intra_O] for rate rule [R5;doublebond_intra_2H_secDe;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C1([CH2])OO[CH]C1=C(23588)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.79781e+09,'s^-1'), n=0.611011, Ea=(94.2067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H;radadd_intra_O] + [R6;doublebond_intra_2H;radadd_intra] for rate rule [R6;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(4)', '[CH2]C(=C)C(=C)C=O(23589)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(21.0272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 20.9 to 21.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C=CO[O](16806)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00267131,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=CO[O](16807)', 'allene(458)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CO[O](16807)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH]=C([CH]O[O])C([CH2])=C(20113)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(CO[O])C([CH2])=C(19489)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C)C(=C)[CH]O[O](23590)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([CH]O[O])C(=C)C(20117)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C(=C)C(=C)CO[O](18231)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C([CH]OO)C([CH2])=C(20114)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([CH2])C(=C)[CH]OO(23591)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=C)C[C]=CO[O](21271)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(T)(20)', '[CH2]C(=C)[C]=CO[O](23592)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['C=C1CCC1=CO[O](22280)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C(=C)C1=COOC1(19219)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C1=COOCC1=C(23593)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(6.42e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSDSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['O(4)', '[CH2]C(=C)C1=COC1(20105)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['O(4)', '[CH2]C1=COCC1=C(23594)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.026e+11,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;C_pri_rad_intra;OO] for rate rule [R4OO_SSD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([CH2])C(=C)[C]O[O](22189)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C](C)C([CH2])=[C]O[O](23595)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]C(=C)[C]([CH2])CO[O](18232)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2][C](C([CH2])=C)C1OO1(23596)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C1([CH2])CC1=CO[O](22341)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(242.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 241.5 to 242.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(=C)C([CH2])=[C]OO(23597)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    products = ['[CH2]C(=C)C(C)=[C]O[O](23598)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C(=[C]O[O])C(=C)C(23599)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=C([CH2])C(C)=CO[O](23600)'],
    products = ['[CH2]C(=C)C(=C)[CH]O[O](20115)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4985',
    isomers = [
        '[CH2]C(=C)C(=C)[CH]O[O](20115)',
    ],
    reactants = [
        ('C=C=CO[O](16806)', 'allene(458)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4985',
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


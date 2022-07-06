species(
    label = '[CH2][CH]CC([CH2])CO[O](742)',
    structure = SMILES('[CH2][CH]CC([CH2])CO[O]'),
    E0 = (470.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,286.686,912.306,2068.72],'cm^-1')),
        HinderedRotor(inertia=(0.0716614,'amu*angstrom^2'), symmetry=1, barrier=(3.38543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716614,'amu*angstrom^2'), symmetry=1, barrier=(3.38543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716614,'amu*angstrom^2'), symmetry=1, barrier=(3.38543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716614,'amu*angstrom^2'), symmetry=1, barrier=(3.38543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716614,'amu*angstrom^2'), symmetry=1, barrier=(3.38543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716614,'amu*angstrom^2'), symmetry=1, barrier=(3.38543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199178,0.0905992,-0.000119451,9.87818e-08,-3.32592e-11,56766.1,37.9845], Tmin=(100,'K'), Tmax=(853.379,'K')), NASAPolynomial(coeffs=[6.77991,0.0469488,-2.0218e-05,3.67731e-09,-2.4706e-13,56109.2,10.0115], Tmin=(853.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'allyl(82)',
    structure = SMILES('[CH2]C=C'),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29615,0.00579223,4.3392e-05,-5.9989e-08,2.33815e-11,18908.2,9.01994], Tmin=(100,'K'), Tmax=(942.181,'K')), NASAPolynomial(coeffs=[8.06863,0.0101837,-2.84795e-06,5.0088e-10,-3.79628e-14,16914.7,-19.5272], Tmin=(942.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'allylperoxy(462)',
    structure = SMILES('C=CCO[O]'),
    E0 = (71.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,262.651],'cm^-1')),
        HinderedRotor(inertia=(0.177302,'amu*angstrom^2'), symmetry=1, barrier=(8.66672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174787,'amu*angstrom^2'), symmetry=1, barrier=(8.66487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3494.07,'J/mol'), sigma=(5.81539,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=545.77 K, Pc=40.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29532,0.0323534,-1.53572e-05,-7.98968e-10,1.98651e-12,8683.79,18.5043], Tmin=(100,'K'), Tmax=(1108.68,'K')), NASAPolynomial(coeffs=[8.87094,0.0173881,-6.95985e-06,1.2773e-09,-8.84653e-14,6687.43,-16.3253], Tmin=(1108.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""allylperoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(CC=C)CO[O](747)',
    structure = SMILES('[CH2]C(CC=C)CO[O]'),
    E0 = (195.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406849,0.0801902,-7.36095e-05,3.92095e-08,-8.8199e-12,23658.6,33.1959], Tmin=(100,'K'), Tmax=(1049.96,'K')), NASAPolynomial(coeffs=[10.873,0.0403179,-1.66471e-05,3.04159e-09,-2.08201e-13,21460.8,-17.8073], Tmin=(1049.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC[CH]CO[O](741)',
    structure = SMILES('[CH2][CH]CC[CH]CO[O]'),
    E0 = (471.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,214.518,810.834,1216.7,1716.14],'cm^-1')),
        HinderedRotor(inertia=(0.124619,'amu*angstrom^2'), symmetry=1, barrier=(3.52154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124619,'amu*angstrom^2'), symmetry=1, barrier=(3.52154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124619,'amu*angstrom^2'), symmetry=1, barrier=(3.52154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124619,'amu*angstrom^2'), symmetry=1, barrier=(3.52154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124619,'amu*angstrom^2'), symmetry=1, barrier=(3.52154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124619,'amu*angstrom^2'), symmetry=1, barrier=(3.52154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.188735,0.0935742,-0.000133684,1.19548e-07,-4.29111e-11,56887.3,38.8229], Tmin=(100,'K'), Tmax=(836.319,'K')), NASAPolynomial(coeffs=[4.98621,0.051442,-2.37048e-05,4.44707e-09,-3.03862e-13,56755.9,20.5469], Tmin=(836.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]CCO[O](1107)',
    structure = SMILES('[CH2][CH]C[CH]CCO[O]'),
    E0 = (465.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418158,0.0920187,-0.00014171,1.36477e-07,-5.09478e-11,56159.3,38.2577], Tmin=(100,'K'), Tmax=(855.636,'K')), NASAPolynomial(coeffs=[1.43654,0.0567431,-2.63744e-05,4.93312e-09,-3.35013e-13,57102,40.0307], Tmin=(855.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])CO[O](744)',
    structure = SMILES('[CH2]C([CH2])C([CH2])CO[O]'),
    E0 = (475.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,210.182,1027.97],'cm^-1')),
        HinderedRotor(inertia=(0.122653,'amu*angstrom^2'), symmetry=1, barrier=(3.36627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122653,'amu*angstrom^2'), symmetry=1, barrier=(3.36627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122653,'amu*angstrom^2'), symmetry=1, barrier=(3.36627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122653,'amu*angstrom^2'), symmetry=1, barrier=(3.36627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122653,'amu*angstrom^2'), symmetry=1, barrier=(3.36627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122653,'amu*angstrom^2'), symmetry=1, barrier=(3.36627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3961.16,'J/mol'), sigma=(7.02744,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=618.72 K, Pc=25.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.109875,0.0919343,-0.000109287,7.62387e-08,-2.14337e-11,57357.1,37.1606], Tmin=(100,'K'), Tmax=(941.11,'K')), NASAPolynomial(coeffs=[11.4776,0.037248,-1.34597e-05,2.21853e-09,-1.40254e-13,55416.8,-16.7597], Tmin=(941.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'O(4)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2][CH]CC([CH2])C[O](1108)',
    structure = SMILES('[CH2][CH]CC([CH2])C[O]'),
    E0 = (473.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,474.928,677.162,2057.93,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0318109,'amu*angstrom^2'), symmetry=1, barrier=(3.62578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318109,'amu*angstrom^2'), symmetry=1, barrier=(3.62578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318109,'amu*angstrom^2'), symmetry=1, barrier=(3.62578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318109,'amu*angstrom^2'), symmetry=1, barrier=(3.62578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318109,'amu*angstrom^2'), symmetry=1, barrier=(3.62578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980149,0.0744452,-9.57014e-05,8.60578e-08,-3.16247e-11,57001,34.5204], Tmin=(100,'K'), Tmax=(845.032,'K')), NASAPolynomial(coeffs=[2.47735,0.0502152,-2.22612e-05,4.11197e-09,-2.78884e-13,57360.1,31.1709], Tmin=(845.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH]CO[O](803)',
    structure = SMILES('[CH2][CH]C[CH]CO[O]'),
    E0 = (495.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,188.467,379.263,3338.83],'cm^-1')),
        HinderedRotor(inertia=(0.029798,'amu*angstrom^2'), symmetry=1, barrier=(0.917953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.029798,'amu*angstrom^2'), symmetry=1, barrier=(0.917953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.029798,'amu*angstrom^2'), symmetry=1, barrier=(0.917953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.029798,'amu*angstrom^2'), symmetry=1, barrier=(0.917953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.029798,'amu*angstrom^2'), symmetry=1, barrier=(0.917953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86742,0.078392,-0.000118294,1.09137e-07,-3.95541e-11,59723.3,34.1422], Tmin=(100,'K'), Tmax=(850.992,'K')), NASAPolynomial(coeffs=[3.99385,0.0431504,-1.99598e-05,3.73107e-09,-2.53575e-13,59935.1,23.9346], Tmin=(850.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2]C([CH2])CO[O](763)',
    structure = SMILES('[CH2]C([CH2])CO[O]'),
    E0 = (320.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.101454,'amu*angstrom^2'), symmetry=1, barrier=(2.33263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0973465,'amu*angstrom^2'), symmetry=1, barrier=(2.23819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096224,'amu*angstrom^2'), symmetry=1, barrier=(2.21238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992394,'amu*angstrom^2'), symmetry=1, barrier=(2.28171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26781,0.0618183,-7.41944e-05,5.33299e-08,-1.54734e-11,38644.1,26.2251], Tmin=(100,'K'), Tmax=(935.864,'K')), NASAPolynomial(coeffs=[8.27545,0.0267561,-9.80543e-06,1.62698e-09,-1.03135e-13,37556.3,-5.92243], Tmin=(935.864,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2][C]CC([CH2])CO[O](1109)',
    structure = SMILES('[CH2][C]CC([CH2])CO[O]'),
    E0 = (724.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230456,0.0894382,-0.000101587,5.9139e-08,-1.02134e-11,87287.1,35.0818], Tmin=(100,'K'), Tmax=(658.817,'K')), NASAPolynomial(coeffs=[10.6033,0.0393441,-1.68677e-05,3.09444e-09,-2.10651e-13,85640.8,-12.7544], Tmin=(658.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(724.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(C[CH][CH2])CO[O](1110)',
    structure = SMILES('[CH]C(C[CH][CH2])CO[O]'),
    E0 = (714.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.166186,0.0915449,-0.000125126,1.04595e-07,-3.56912e-11,86009.2,37.3249], Tmin=(100,'K'), Tmax=(824.177,'K')), NASAPolynomial(coeffs=[7.55627,0.0448085,-2.02823e-05,3.78594e-09,-2.58589e-13,85160.3,5.34074], Tmin=(824.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC([CH2])CO[O](1111)',
    structure = SMILES('[CH][CH]CC([CH2])CO[O]'),
    E0 = (713.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0931629,0.0931155,-0.000129187,1.07033e-07,-3.55147e-11,85992.2,37.5823], Tmin=(100,'K'), Tmax=(864.195,'K')), NASAPolynomial(coeffs=[7.95923,0.0430318,-1.85198e-05,3.3509e-09,-2.2383e-13,85143.3,3.73563], Tmin=(864.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC(CO[O])C1(1091)',
    structure = SMILES('[CH2]C1CC(CO[O])C1'),
    E0 = (216.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.746495,0.0586114,-4.6078e-06,-3.33231e-08,1.64635e-11,26226.9,30.0911], Tmin=(100,'K'), Tmax=(989.525,'K')), NASAPolynomial(coeffs=[14.0724,0.034315,-1.26045e-05,2.26548e-09,-1.58031e-13,22141.9,-41.374], Tmin=(989.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC1COOC1(1112)',
    structure = SMILES('[CH2][CH]CC1COOC1'),
    E0 = (194.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931753,0.0568895,-4.1201e-07,-4.12238e-08,2.16933e-11,23560.2,30.2491], Tmin=(100,'K'), Tmax=(899.636,'K')), NASAPolynomial(coeffs=[12.6048,0.0338324,-1.00606e-05,1.56478e-09,-1.00812e-13,20292.7,-31.3191], Tmin=(899.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1COOC([CH2])C1(1113)',
    structure = SMILES('[CH2]C1COOC([CH2])C1'),
    E0 = (193.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.733385,0.052145,3.7122e-05,-9.79837e-08,4.67538e-11,23390.2,27.2857], Tmin=(100,'K'), Tmax=(889.1,'K')), NASAPolynomial(coeffs=[19.4331,0.0219102,-2.79312e-06,1.22638e-10,-2.35957e-15,17934.9,-72.7107], Tmin=(889.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(Isobutyl) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][CH]CC1COC1(1114)',
    structure = SMILES('[CH2][CH]CC1COC1'),
    E0 = (222.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52787,0.0470149,5.24765e-06,-3.75728e-08,1.8689e-11,26915.6,26.6315], Tmin=(100,'K'), Tmax=(881.458,'K')), NASAPolynomial(coeffs=[8.65037,0.0356496,-1.10731e-05,1.7425e-09,-1.11308e-13,24845.8,-11.4497], Tmin=(881.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1COC([CH2])C1(1115)',
    structure = SMILES('[CH2]C1COC([CH2])C1'),
    E0 = (147.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25817,0.0470312,2.39941e-05,-7.18964e-08,3.51739e-11,17885.9,22.0696], Tmin=(100,'K'), Tmax=(870.697,'K')), NASAPolynomial(coeffs=[14.1441,0.0258681,-5.07224e-06,5.29537e-10,-2.68082e-14,14200.2,-46.5924], Tmin=(870.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(Isobutyl) + radical(CJC(C)OC)"""),
)

species(
    label = 'HO2(9)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26737e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.914,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11273e-07,6.34955e-11,-4.86383e-15,83.4204,3.0934], Tmin=(923.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2][CH]CC([CH2])=C(643)',
    structure = SMILES('[CH2][CH]CC([CH2])=C'),
    E0 = (471.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,246.856,1775.46],'cm^-1')),
        HinderedRotor(inertia=(0.0051858,'amu*angstrom^2'), symmetry=1, barrier=(11.5977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.34837e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701505,'amu*angstrom^2'), symmetry=1, barrier=(30.3593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01076,'amu*angstrom^2'), symmetry=1, barrier=(87.6001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42509,0.0531145,-3.31274e-05,1.07407e-08,-1.45659e-12,56843.8,26.4331], Tmin=(100,'K'), Tmax=(1641.34,'K')), NASAPolynomial(coeffs=[10.9943,0.029794,-1.1815e-05,2.08414e-09,-1.38063e-13,53702.6,-24.4744], Tmin=(1641.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC(=C)CO[O](1116)',
    structure = SMILES('[CH2]CCC(=C)CO[O]'),
    E0 = (196.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615008,0.0790237,-7.33925e-05,4.09157e-08,-9.94879e-12,23723,32.307], Tmin=(100,'K'), Tmax=(958.034,'K')), NASAPolynomial(coeffs=[8.77176,0.044968,-2.0072e-05,3.81215e-09,-2.66703e-13,22160.1,-6.69508], Tmin=(958.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC(C)CO[O](1117)',
    structure = SMILES('[CH2]C=CC(C)CO[O]'),
    E0 = (130.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.307703,0.0775196,-6.13089e-05,2.58959e-08,-4.52637e-12,15877.2,31.6407], Tmin=(100,'K'), Tmax=(1334.85,'K')), NASAPolynomial(coeffs=[14.0896,0.0362212,-1.49011e-05,2.71848e-09,-1.85574e-13,12197.9,-38.8294], Tmin=(1334.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC(=C)COO(1118)',
    structure = SMILES('[CH2][CH]CC(=C)COO'),
    E0 = (238.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.458379,0.0839071,-9.20333e-05,6.56328e-08,-2.05779e-11,28831.6,34.8455], Tmin=(100,'K'), Tmax=(756.314,'K')), NASAPolynomial(coeffs=[7.09098,0.0488258,-2.24508e-05,4.29311e-09,-3.00354e-13,27828.4,4.70002], Tmin=(756.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])COO(1119)',
    structure = SMILES('[CH2]C=CC([CH2])COO'),
    E0 = (183.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0294339,0.0838602,-7.27985e-05,3.38794e-08,-6.4335e-12,22274.2,33.8745], Tmin=(100,'K'), Tmax=(1252.97,'K')), NASAPolynomial(coeffs=[15.412,0.0345644,-1.37835e-05,2.47925e-09,-1.68314e-13,18404.7,-44.1037], Tmin=(1252.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC(=C)CO[O](1120)',
    structure = SMILES('[CH2][CH]CC(=C)CO[O]'),
    E0 = (390.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,221.487,821.487,1716.24],'cm^-1')),
        HinderedRotor(inertia=(0.125153,'amu*angstrom^2'), symmetry=1, barrier=(3.529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125153,'amu*angstrom^2'), symmetry=1, barrier=(3.529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125153,'amu*angstrom^2'), symmetry=1, barrier=(3.529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125153,'amu*angstrom^2'), symmetry=1, barrier=(3.529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125153,'amu*angstrom^2'), symmetry=1, barrier=(3.529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543062,0.0826513,-0.000106006,8.86231e-08,-3.0958e-11,47109.6,35.286], Tmin=(100,'K'), Tmax=(808.987,'K')), NASAPolynomial(coeffs=[5.94969,0.0463539,-2.09697e-05,3.93271e-09,-2.70051e-13,46547.8,12.2826], Tmin=(808.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])CO[O](1121)',
    structure = SMILES('[CH2]C=CC([CH2])CO[O]'),
    E0 = (335.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,199.814,4000],'cm^-1')),
        HinderedRotor(inertia=(0.36512,'amu*angstrom^2'), symmetry=1, barrier=(10.3502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.13979,'amu*angstrom^2'), symmetry=1, barrier=(117.356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0042215,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.13947,'amu*angstrom^2'), symmetry=1, barrier=(117.355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.13835,'amu*angstrom^2'), symmetry=1, barrier=(117.35,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341696,0.079008,-7.31871e-05,3.78478e-08,-8.07076e-12,40540.1,33.3023], Tmin=(100,'K'), Tmax=(1118.75,'K')), NASAPolynomial(coeffs=[12.7215,0.0347457,-1.38419e-05,2.48434e-09,-1.68448e-13,37770.1,-27.8126], Tmin=(1118.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH2](497)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.000863049,'amu*angstrom^2'), symmetry=1, barrier=(2.40754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000864751,'amu*angstrom^2'), symmetry=1, barrier=(2.41365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42318,0.0135015,1.25385e-06,-4.52647e-09,1.21417e-12,58330.2,15.4092], Tmin=(100,'K'), Tmax=(1578.37,'K')), NASAPolynomial(coeffs=[5.16145,0.0152713,-6.29635e-06,1.14122e-09,-7.61383e-14,57012.3,3.79304], Tmin=(1578.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]O[O](46)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (205.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2028.66,2028.67],'cm^-1')),
        HinderedRotor(inertia=(0.221535,'amu*angstrom^2'), symmetry=1, barrier=(13.5753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43223,0.0119348,-9.20111e-06,4.32805e-09,-8.54519e-13,24777.3,10.782], Tmin=(100,'K'), Tmax=(1197.14,'K')), NASAPolynomial(coeffs=[5.14053,0.00622677,-2.04895e-06,3.45073e-10,-2.27385e-14,24368.2,2.23311], Tmin=(1197.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = '[CH2][CH]CC=C(508)',
    structure = SMILES('[CH2][CH]CC=C'),
    E0 = (359.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1243.57],'cm^-1')),
        HinderedRotor(inertia=(0.00331557,'amu*angstrom^2'), symmetry=1, barrier=(3.65555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160798,'amu*angstrom^2'), symmetry=1, barrier=(3.69707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161411,'amu*angstrom^2'), symmetry=1, barrier=(3.71116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13435,0.0369817,-1.51234e-05,1.33725e-09,3.77008e-13,43295.2,23.5599], Tmin=(100,'K'), Tmax=(1532.36,'K')), NASAPolynomial(coeffs=[8.68459,0.025936,-1.02358e-05,1.78847e-09,-1.17139e-13,40577.1,-13.1557], Tmin=(1532.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CO[O](461)',
    structure = SMILES('[CH2][CH]CO[O]'),
    E0 = (356.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.112977,'amu*angstrom^2'), symmetry=1, barrier=(2.59756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112908,'amu*angstrom^2'), symmetry=1, barrier=(2.59597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113038,'amu*angstrom^2'), symmetry=1, barrier=(2.59896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19834,0.0429895,-5.67485e-05,4.73876e-08,-1.62792e-11,42972,23.4058], Tmin=(100,'K'), Tmax=(810.767,'K')), NASAPolynomial(coeffs=[5.55675,0.0219523,-9.56083e-06,1.78939e-09,-1.23027e-13,42574.3,8.8136], Tmin=(810.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCJCOOH) + radical(RCCJ)"""),
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
    label = '[CH2][CH]CC([CH2])[CH2](484)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH2]'),
    E0 = (612.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1170.47,1171.2],'cm^-1')),
        HinderedRotor(inertia=(0.169767,'amu*angstrom^2'), symmetry=1, barrier=(3.90329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00400723,'amu*angstrom^2'), symmetry=1, barrier=(3.90322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169659,'amu*angstrom^2'), symmetry=1, barrier=(3.90079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169917,'amu*angstrom^2'), symmetry=1, barrier=(3.90674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00401056,'amu*angstrom^2'), symmetry=1, barrier=(3.90297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64314,0.0551022,-3.22884e-05,-1.61575e-09,1.00735e-11,73724.5,29.7062], Tmin=(100,'K'), Tmax=(647.823,'K')), NASAPolynomial(coeffs=[5.92384,0.0380634,-1.45841e-05,2.54566e-09,-1.69302e-13,72972.8,9.39149], Tmin=(647.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH][CH2](510)',
    structure = SMILES('[CH2][CH]C[CH][CH2]'),
    E0 = (631.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1407.99,1417.06],'cm^-1')),
        HinderedRotor(inertia=(0.00560336,'amu*angstrom^2'), symmetry=1, barrier=(7.95035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00563889,'amu*angstrom^2'), symmetry=1, barrier=(7.95031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00557143,'amu*angstrom^2'), symmetry=1, barrier=(7.9482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00562946,'amu*angstrom^2'), symmetry=1, barrier=(7.93664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26906,0.0435841,-4.87505e-05,4.66768e-08,-1.86942e-11,75984.9,26.4118], Tmin=(100,'K'), Tmax=(836.615,'K')), NASAPolynomial(coeffs=[0.578549,0.0388266,-1.7199e-05,3.18971e-09,-2.17339e-13,76717.2,36.9515], Tmin=(836.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])CO[O](1122)',
    structure = SMILES('[CH2][CH]C[C]([CH2])CO[O]'),
    E0 = (661.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,375.864,3120.76,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0217342,'amu*angstrom^2'), symmetry=1, barrier=(3.50732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217342,'amu*angstrom^2'), symmetry=1, barrier=(3.50732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217342,'amu*angstrom^2'), symmetry=1, barrier=(3.50732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217342,'amu*angstrom^2'), symmetry=1, barrier=(3.50732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217342,'amu*angstrom^2'), symmetry=1, barrier=(3.50732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217342,'amu*angstrom^2'), symmetry=1, barrier=(3.50732,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782773,0.0773432,-9.68902e-05,8.18106e-08,-2.86469e-11,79706.8,40.7342], Tmin=(100,'K'), Tmax=(840.659,'K')), NASAPolynomial(coeffs=[4.61752,0.0467901,-2.04146e-05,3.749e-09,-2.53785e-13,79496.9,25.4858], Tmin=(840.659,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])CO[O](1123)',
    structure = SMILES('[CH2][CH][CH]C([CH2])CO[O]'),
    E0 = (665.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,225.173,1963.88,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0716693,'amu*angstrom^2'), symmetry=1, barrier=(2.69402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716693,'amu*angstrom^2'), symmetry=1, barrier=(2.69402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716693,'amu*angstrom^2'), symmetry=1, barrier=(2.69402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716693,'amu*angstrom^2'), symmetry=1, barrier=(2.69402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716693,'amu*angstrom^2'), symmetry=1, barrier=(2.69402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716693,'amu*angstrom^2'), symmetry=1, barrier=(2.69402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.38184,0.0895015,-0.000131407,1.17961e-07,-4.17231e-11,80154.8,39.8836], Tmin=(100,'K'), Tmax=(862.341,'K')), NASAPolynomial(coeffs=[4.7245,0.0478081,-2.13985e-05,3.93546e-09,-2.6483e-13,80207,24.2218], Tmin=(862.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])[CH]O[O](1124)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH]O[O]'),
    E0 = (659.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,304.889,934.68,2144.55],'cm^-1')),
        HinderedRotor(inertia=(0.063266,'amu*angstrom^2'), symmetry=1, barrier=(3.38041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.063266,'amu*angstrom^2'), symmetry=1, barrier=(3.38041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.063266,'amu*angstrom^2'), symmetry=1, barrier=(3.38041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.063266,'amu*angstrom^2'), symmetry=1, barrier=(3.38041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.063266,'amu*angstrom^2'), symmetry=1, barrier=(3.38041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.063266,'amu*angstrom^2'), symmetry=1, barrier=(3.38041,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.170836,0.0937268,-0.000137683,1.20619e-07,-4.15587e-11,79446,38.2347], Tmin=(100,'K'), Tmax=(871.776,'K')), NASAPolynomial(coeffs=[6.06399,0.0460751,-2.02263e-05,3.6766e-09,-2.45411e-13,79201.8,15.1045], Tmin=(871.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJC) + radical(CCsJOOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C](C)CO[O](1125)',
    structure = SMILES('[CH2][CH]C[C](C)CO[O]'),
    E0 = (456.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.875045,0.0743361,-7.97815e-05,6.36627e-08,-2.29397e-11,55038.7,38.624], Tmin=(100,'K'), Tmax=(764.414,'K')), NASAPolynomial(coeffs=[4.16353,0.0513199,-2.32192e-05,4.39295e-09,-3.04756e-13,54705.7,24.7524], Tmin=(764.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C([CH2])CO[O](1126)',
    structure = SMILES('[CH2]C[CH]C([CH2])CO[O]'),
    E0 = (470.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,244.763,883.931,1958.11],'cm^-1')),
        HinderedRotor(inertia=(0.100811,'amu*angstrom^2'), symmetry=1, barrier=(3.4715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100811,'amu*angstrom^2'), symmetry=1, barrier=(3.4715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100811,'amu*angstrom^2'), symmetry=1, barrier=(3.4715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100811,'amu*angstrom^2'), symmetry=1, barrier=(3.4715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100811,'amu*angstrom^2'), symmetry=1, barrier=(3.4715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100811,'amu*angstrom^2'), symmetry=1, barrier=(3.4715,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.224996,0.0887445,-0.000109777,8.61987e-08,-2.84638e-11,56778,37.7149], Tmin=(100,'K'), Tmax=(802.23,'K')), NASAPolynomial(coeffs=[7.7344,0.0461,-2.03144e-05,3.77088e-09,-2.57833e-13,55740.5,4.18442], Tmin=(802.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])[CH]OO(1127)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH]OO'),
    E0 = (507.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0628138,0.0969478,-0.000131646,1.09615e-07,-3.71367e-11,61174.3,38.3155], Tmin=(100,'K'), Tmax=(837.308,'K')), NASAPolynomial(coeffs=[7.61091,0.0478098,-2.12618e-05,3.9281e-09,-2.6645e-13,60326.7,5.26898], Tmin=(837.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C(C)CO[O](1128)',
    structure = SMILES('[CH2][CH][CH]C(C)CO[O]'),
    E0 = (460.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.42536,0.0871155,-0.000116676,1.03146e-07,-3.75164e-11,55488.8,37.9452], Tmin=(100,'K'), Tmax=(817.482,'K')), NASAPolynomial(coeffs=[4.45293,0.05201,-2.40065e-05,4.53159e-09,-3.11749e-13,55344.8,22.4731], Tmin=(817.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC(C)[CH]O[O](1129)',
    structure = SMILES('[CH2][CH]CC(C)[CH]O[O]'),
    E0 = (454.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216192,0.0913107,-0.000122794,1.05493e-07,-3.71525e-11,54780,36.2901], Tmin=(100,'K'), Tmax=(825.363,'K')), NASAPolynomial(coeffs=[5.81368,0.0502401,-2.28126e-05,4.26755e-09,-2.91897e-13,54330.9,13.2369], Tmin=(825.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJC) + radical(CCsJOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[C]([CH2])CO[O](1130)',
    structure = SMILES('[CH2]CC[C]([CH2])CO[O]'),
    E0 = (467.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.781149,0.0746135,-6.76108e-05,3.88452e-08,-9.92365e-12,56323.4,38.0175], Tmin=(100,'K'), Tmax=(909.99,'K')), NASAPolynomial(coeffs=[7.41382,0.0454587,-1.95532e-05,3.63795e-09,-2.51286e-13,55116.2,6.64422], Tmin=(909.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH][CH]C)CO[O](886)',
    structure = SMILES('[CH2]C([CH][CH]C)CO[O]'),
    E0 = (460.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,374.458,515.766,3185.72],'cm^-1')),
        HinderedRotor(inertia=(0.0458994,'amu*angstrom^2'), symmetry=1, barrier=(1.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458994,'amu*angstrom^2'), symmetry=1, barrier=(1.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458994,'amu*angstrom^2'), symmetry=1, barrier=(1.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458994,'amu*angstrom^2'), symmetry=1, barrier=(1.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458994,'amu*angstrom^2'), symmetry=1, barrier=(1.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0458994,'amu*angstrom^2'), symmetry=1, barrier=(1.8359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.360753,0.0885717,-0.000120252,1.04798e-07,-3.69132e-11,55471.4,38.1733], Tmin=(100,'K'), Tmax=(849.464,'K')), NASAPolynomial(coeffs=[4.85911,0.0502285,-2.22415e-05,4.09601e-09,-2.76947e-13,55326.3,20.8497], Tmin=(849.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCC([CH2])[CH]O[O](1131)',
    structure = SMILES('[CH2]CCC([CH2])[CH]O[O]'),
    E0 = (465.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00912308,0.0930275,-0.000116246,8.90704e-08,-2.83586e-11,56069.5,36.0833], Tmin=(100,'K'), Tmax=(814.908,'K')), NASAPolynomial(coeffs=[9.11135,0.0443006,-1.91028e-05,3.5025e-09,-2.37613e-13,54720.4,-5.14192], Tmin=(814.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C[CH]C)CO[O](885)',
    structure = SMILES('[CH2][C](C[CH]C)CO[O]'),
    E0 = (456.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,242.365,2061.85,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0506566,'amu*angstrom^2'), symmetry=1, barrier=(2.23218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506566,'amu*angstrom^2'), symmetry=1, barrier=(2.23218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506566,'amu*angstrom^2'), symmetry=1, barrier=(2.23218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506566,'amu*angstrom^2'), symmetry=1, barrier=(2.23218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506566,'amu*angstrom^2'), symmetry=1, barrier=(2.23218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506566,'amu*angstrom^2'), symmetry=1, barrier=(2.23218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03409,0.0549231,3.67829e-05,-2.16027e-07,2.07767e-10,54970.1,34.7935], Tmin=(100,'K'), Tmax=(417.724,'K')), NASAPolynomial(coeffs=[5.01057,0.0487098,-2.09412e-05,3.82978e-09,-2.58958e-13,54526.9,20.7044], Tmin=(417.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])COO(1132)',
    structure = SMILES('[CH2][CH]C[C]([CH2])COO'),
    E0 = (509.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,354.675,1245.77,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0395711,'amu*angstrom^2'), symmetry=1, barrier=(1.55937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36012,0.0679718,-2.67327e-05,-5.92022e-08,6.68434e-11,61400.4,38.0675], Tmin=(100,'K'), Tmax=(495.936,'K')), NASAPolynomial(coeffs=[5.82053,0.0491779,-2.1857e-05,4.10214e-09,-2.8361e-13,60746.7,17.5465], Tmin=(495.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CJCOOH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]O[O])C[CH]C(887)',
    structure = SMILES('[CH2]C([CH]O[O])C[CH]C'),
    E0 = (454.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.1482,0.0928151,-0.000126587,1.07521e-07,-3.67669e-11,54762.7,36.5299], Tmin=(100,'K'), Tmax=(859.494,'K')), NASAPolynomial(coeffs=[6.21009,0.0484754,-2.10574e-05,3.8343e-09,-2.5729e-13,54316.4,11.6682], Tmin=(859.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJC) + radical(CCsJOOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])COO(1133)',
    structure = SMILES('[CH2][CH][CH]C([CH2])COO'),
    E0 = (513.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,242.874,880.388,1942.99],'cm^-1')),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0975941,'amu*angstrom^2'), symmetry=1, barrier=(3.30902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.149024,0.0927163,-0.000125375,1.07022e-07,-3.73673e-11,61883,39.9613], Tmin=(100,'K'), Tmax=(827.632,'K')), NASAPolynomial(coeffs=[6.24987,0.0495805,-2.24562e-05,4.19228e-09,-2.86315e-13,61340.6,14.5068], Tmin=(827.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (470.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (470.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (630.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (630.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (632.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (716.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (933.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (933.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (936.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (925.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (925.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (479.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (478.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (479.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (553.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (525.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (610.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (534.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (534.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (495.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (495.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (602.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (555.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (570.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (601.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (550.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (603.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (841.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (837.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (873.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (877.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (871.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (573.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (607.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (636.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (588.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (588.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (587.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (612.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (561.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (602.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (593.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (525.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (555.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['allyl(82)', 'allylperoxy(462)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C(CC=C)CO[O](747)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]CC[CH]CO[O](741)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]C[CH]CCO[O](1107)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2][CH]CC([CH2])C[O](1108)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][CH]C[CH]CO[O](803)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH2](502)', '[CH2]C([CH2])CO[O](763)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]CC([CH2])CO[O](1109)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C(C[CH][CH2])CO[O](1110)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH]CC([CH2])CO[O](1111)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C1CC(CO[O])C1(1091)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]CC1COOC1(1112)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C1COOC([CH2])C1(1113)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['O(4)', '[CH2][CH]CC1COC1(1114)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['O(4)', '[CH2]C1COC([CH2])C1(1115)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO_SSS;Y_rad_intra;OO] for rate rule [R4OO_SSS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['HO2(9)', '[CH2][CH]CC([CH2])=C(643)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]CCC(=C)CO[O](1116)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C=CC(C)CO[O](1117)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]CC(=C)COO(1118)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C=CC([CH2])COO(1119)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH]CC(=C)CO[O](1120)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C=CC([CH2])CO[O](1121)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][CH2](497)', 'allylperoxy(462)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.08019,'m^3/(mol*s)'), n=2.05449, Ea=(13.6169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;CJ] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]O[O](46)', '[CH2][CH]CC=C(508)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00140605,'m^3/(mol*s)'), n=2.41796, Ea=(36.1629,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CsJ] for rate rule [Cds-CsH_Cds-HH;CsJ-OsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['allyl(82)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O2(2)', '[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.4339e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(616.814,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]O[O](46)', '[CH2][CH]C[CH][CH2](510)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][CH]C[C]([CH2])CO[O](1122)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][CH][CH]C([CH2])CO[O](1123)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH]CC([CH2])[CH]O[O](1124)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]C[C](C)CO[O](1125)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C[CH]C([CH2])CO[O](1126)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]CC([CH2])[CH]OO(1127)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 242 used for R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH][CH]C(C)CO[O](1128)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][CH]CC(C)[CH]O[O](1129)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]CC[C]([CH2])CO[O](1130)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C([CH][CH]C)CO[O](886)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CCC([CH2])[CH]O[O](1131)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(55800,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2][C](C[CH]C)CO[O](885)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C[C]([CH2])COO(1132)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]CC([CH2])CO[O](742)'],
    products = ['[CH2]C([CH]O[O])C[CH]C(887)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(126713,'s^-1'), n=1.75034, Ea=(55.0506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH][CH]C([CH2])COO(1133)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #312',
    isomers = [
        '[CH2][CH]CC([CH2])CO[O](742)',
    ],
    reactants = [
        ('allyl(82)', 'allylperoxy(462)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #312',
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


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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.109875,0.0919343,-0.000109287,7.62387e-08,-2.14337e-11,57357.1,37.1606], Tmin=(100,'K'), Tmax=(941.11,'K')), NASAPolynomial(coeffs=[11.4776,0.037248,-1.34597e-05,2.21853e-09,-1.40254e-13,55416.8,-16.7597], Tmin=(941.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH2])C[CH]CO[O](743)',
    structure = SMILES('[CH2]C([CH2])C[CH]CO[O]'),
    E0 = (476.703,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0603537,0.094132,-0.000120429,9.23531e-08,-2.87563e-11,57475.7,37.7883], Tmin=(100,'K'), Tmax=(872.949,'K')), NASAPolynomial(coeffs=[9.62546,0.0418482,-1.70113e-05,3.00417e-09,-1.98407e-13,56085.7,-5.89963], Tmin=(872.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH]CCO[O](1061)',
    structure = SMILES('[CH2]C([CH2])[CH]CCO[O]'),
    E0 = (470.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.129647,0.0906022,-0.000114958,9.02608e-08,-2.90574e-11,56761.9,37.3577], Tmin=(100,'K'), Tmax=(856.44,'K')), NASAPolynomial(coeffs=[8.20444,0.0442029,-1.84799e-05,3.31832e-09,-2.21591e-13,55697.3,1.51264], Tmin=(856.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

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
    label = '[CH2]C([CH2])C([CH2])C[O](1062)',
    structure = SMILES('[CH2]C([CH2])C([CH2])C[O]'),
    E0 = (477.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,361,482.261,3249.12],'cm^-1')),
        HinderedRotor(inertia=(0.0250514,'amu*angstrom^2'), symmetry=1, barrier=(2.76953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0250514,'amu*angstrom^2'), symmetry=1, barrier=(2.76953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0250514,'amu*angstrom^2'), symmetry=1, barrier=(2.76953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0250514,'amu*angstrom^2'), symmetry=1, barrier=(2.76953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0250514,'amu*angstrom^2'), symmetry=1, barrier=(2.76953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.711564,0.075259,-8.34835e-05,6.04582e-08,-1.82879e-11,57590.2,33.5541], Tmin=(100,'K'), Tmax=(903.217,'K')), NASAPolynomial(coeffs=[7.12779,0.0406,-1.55545e-05,2.66577e-09,-1.73146e-13,56685.9,4.66259], Tmin=(903.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH2])[CH]CO[O](1063)',
    structure = SMILES('[CH2]C([CH2])[CH]CO[O]'),
    E0 = (500.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1137.81],'cm^-1')),
        HinderedRotor(inertia=(0.116023,'amu*angstrom^2'), symmetry=1, barrier=(2.66759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114664,'amu*angstrom^2'), symmetry=1, barrier=(2.63635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114794,'amu*angstrom^2'), symmetry=1, barrier=(2.63934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116491,'amu*angstrom^2'), symmetry=1, barrier=(2.67836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115728,'amu*angstrom^2'), symmetry=1, barrier=(2.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612123,0.0790389,-0.000105444,8.26513e-08,-2.58103e-11,60311.9,33.129], Tmin=(100,'K'), Tmax=(902.697,'K')), NASAPolynomial(coeffs=[8.60185,0.0336101,-1.32974e-05,2.29555e-09,-1.48736e-13,59277.9,-2.33655], Tmin=(902.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])CO[O](778)',
    structure = SMILES('[CH2][CH]C([CH2])CO[O]'),
    E0 = (494.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1318.49,1318.67],'cm^-1')),
        HinderedRotor(inertia=(0.154844,'amu*angstrom^2'), symmetry=1, barrier=(3.56017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154875,'amu*angstrom^2'), symmetry=1, barrier=(3.56088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154792,'amu*angstrom^2'), symmetry=1, barrier=(3.55897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154918,'amu*angstrom^2'), symmetry=1, barrier=(3.56187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154747,'amu*angstrom^2'), symmetry=1, barrier=(3.55794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.878934,0.0738813,-9.56361e-05,7.76087e-08,-2.59767e-11,59615,33.1212], Tmin=(100,'K'), Tmax=(830.661,'K')), NASAPolynomial(coeffs=[6.81222,0.0376822,-1.64938e-05,3.03648e-09,-2.05988e-13,58892.4,7.1814], Tmin=(830.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH]C(CO[O])C([CH2])[CH2](1064)',
    structure = SMILES('[CH]C(CO[O])C([CH2])[CH2]'),
    E0 = (718.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140706,0.089041,-9.87797e-05,5.60234e-08,-9.85625e-12,86588.1,35.513], Tmin=(100,'K'), Tmax=(710.738,'K')), NASAPolynomial(coeffs=[11.9663,0.0356441,-1.38535e-05,2.40862e-09,-1.58768e-13,84574.8,-19.8384], Tmin=(710.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(718.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C([CH2])CO[O](1065)',
    structure = SMILES('[CH]C([CH2])C([CH2])CO[O]'),
    E0 = (718.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140706,0.089041,-9.87797e-05,5.60234e-08,-9.85625e-12,86588.1,36.2061], Tmin=(100,'K'), Tmax=(710.738,'K')), NASAPolynomial(coeffs=[11.9663,0.0356441,-1.38535e-05,2.40862e-09,-1.58768e-13,84574.8,-19.1453], Tmin=(710.738,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(718.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(CO[O])C1CC1(1066)',
    structure = SMILES('[CH2]C(CO[O])C1CC1'),
    E0 = (221.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443653,0.0678156,-3.15503e-05,-6.63406e-09,7.48389e-12,26767.1,31.8568], Tmin=(100,'K'), Tmax=(999.956,'K')), NASAPolynomial(coeffs=[14.5633,0.0330839,-1.20762e-05,2.13406e-09,-1.46326e-13,22855.9,-41.699], Tmin=(999.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CCC1CO[O](1067)',
    structure = SMILES('[CH2]C1CCC1CO[O]'),
    E0 = (216.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.746495,0.0586114,-4.6078e-06,-3.33231e-08,1.64635e-11,26226.9,30.0911], Tmin=(100,'K'), Tmax=(989.525,'K')), NASAPolynomial(coeffs=[14.0724,0.034315,-1.26045e-05,2.26548e-09,-1.58031e-13,22141.9,-41.374], Tmin=(989.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C1COOC1(1068)',
    structure = SMILES('[CH2]C([CH2])C1COOC1'),
    E0 = (196.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771287,0.0563342,1.70554e-05,-7.43328e-08,3.85664e-11,23742.3,28.901], Tmin=(100,'K'), Tmax=(862.357,'K')), NASAPolynomial(coeffs=[16.9659,0.0247293,-3.65767e-06,1.91845e-10,-1.24401e-15,19331.2,-56.2111], Tmin=(862.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1COOCC1[CH2](1069)',
    structure = SMILES('[CH2]C1COOCC1[CH2]'),
    E0 = (188.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891694,0.0499949,4.00134e-05,-1.00505e-07,4.82893e-11,22748.5,25.8794], Tmin=(100,'K'), Tmax=(872.852,'K')), NASAPolynomial(coeffs=[18.1524,0.0227732,-2.35926e-06,-4.86135e-11,1.35325e-14,17759.1,-66.3668], Tmin=(872.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C1COC1(1070)',
    structure = SMILES('[CH2]C([CH2])C1COC1'),
    E0 = (224.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36667,0.0464465,2.28924e-05,-7.11716e-08,3.59262e-11,27097.7,25.2872], Tmin=(100,'K'), Tmax=(844.084,'K')), NASAPolynomial(coeffs=[13.0788,0.0264318,-4.60396e-06,3.53899e-10,-1.044e-14,23856.3,-36.7204], Tmin=(844.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1COCC1[CH2](1071)',
    structure = SMILES('[CH2]C1COCC1[CH2]'),
    E0 = (145.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68825,0.0350305,5.6435e-05,-1.05364e-07,4.75916e-11,17653.2,21.9295], Tmin=(100,'K'), Tmax=(865.511,'K')), NASAPolynomial(coeffs=[13.3584,0.0254391,-3.79168e-06,2.19408e-10,-3.90152e-15,13972.2,-42.2812], Tmin=(865.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C(=C)C([CH2])[CH2](635)',
    structure = SMILES('[CH2]C(=C)C([CH2])[CH2]'),
    E0 = (478.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,1525.66],'cm^-1')),
        HinderedRotor(inertia=(0.0127106,'amu*angstrom^2'), symmetry=1, barrier=(8.96111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144501,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.865711,'amu*angstrom^2'), symmetry=1, barrier=(71.4941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862351,'amu*angstrom^2'), symmetry=1, barrier=(71.4907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20734,0.0531304,-1.82483e-05,-1.71517e-08,1.2075e-11,57660.8,25.3167], Tmin=(100,'K'), Tmax=(915.388,'K')), NASAPolynomial(coeffs=[12.7289,0.0246933,-7.55107e-06,1.20412e-09,-7.88685e-14,54633.5,-34.2634], Tmin=(915.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(CO[O])C(=C)C(829)',
    structure = SMILES('[CH2]C(CO[O])C(=C)C'),
    E0 = (185.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.159759,'amu*angstrom^2'), symmetry=1, barrier=(3.67318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159701,'amu*angstrom^2'), symmetry=1, barrier=(3.67185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159758,'amu*angstrom^2'), symmetry=1, barrier=(3.67315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159743,'amu*angstrom^2'), symmetry=1, barrier=(3.6728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.830307,'amu*angstrom^2'), symmetry=1, barrier=(19.0904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283534,0.0841912,-8.52385e-05,5.09803e-08,-1.2809e-11,22456.6,32.1889], Tmin=(100,'K'), Tmax=(950.161,'K')), NASAPolynomial(coeffs=[10.5501,0.0409701,-1.7005e-05,3.1044e-09,-2.12033e-13,20505.7,-16.8163], Tmin=(950.161,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C(=C)CO[O](828)',
    structure = SMILES('[CH2]C(C)C(=C)CO[O]'),
    E0 = (192.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,264.505,2569.54],'cm^-1')),
        HinderedRotor(inertia=(0.188653,'amu*angstrom^2'), symmetry=1, barrier=(9.36605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00240954,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188655,'amu*angstrom^2'), symmetry=1, barrier=(9.36601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188653,'amu*angstrom^2'), symmetry=1, barrier=(9.36603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45693,'amu*angstrom^2'), symmetry=1, barrier=(72.3319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41898,0.0795436,-7.27268e-05,3.87507e-08,-8.71511e-12,23258,32.7498], Tmin=(100,'K'), Tmax=(1051.17,'K')), NASAPolynomial(coeffs=[10.8147,0.0399855,-1.62791e-05,2.95141e-09,-2.0109e-13,21072.4,-17.9225], Tmin=(1051.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C(C)CO[O](1072)',
    structure = SMILES('[CH2]C(=C)C(C)CO[O]'),
    E0 = (132.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.173765,0.0812041,-6.8751e-05,3.13593e-08,-5.89659e-12,16020.4,30.3792], Tmin=(100,'K'), Tmax=(1253.63,'K')), NASAPolynomial(coeffs=[14.2754,0.0362091,-1.49128e-05,2.72844e-09,-1.86946e-13,12484.8,-40.8406], Tmin=(1253.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)COO(1073)',
    structure = SMILES('[CH2]C([CH2])C(=C)COO'),
    E0 = (245.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0626392,0.0861348,-8.51548e-05,4.79576e-08,-1.11211e-11,29655.7,34.3571], Tmin=(100,'K'), Tmax=(1035.91,'K')), NASAPolynomial(coeffs=[12.6682,0.0374603,-1.46741e-05,2.59927e-09,-1.74599e-13,27044,-26.9023], Tmin=(1035.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])COO(1074)',
    structure = SMILES('[CH2]C(=C)C([CH2])COO'),
    E0 = (185.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.172821,0.0876556,-8.06118e-05,3.9781e-08,-7.96719e-12,22417.8,32.6468], Tmin=(100,'K'), Tmax=(1196.49,'K')), NASAPolynomial(coeffs=[15.8142,0.0342088,-1.36066e-05,2.44625e-09,-1.6622e-13,18592.2,-47.3491], Tmin=(1196.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)CO[O](1075)',
    structure = SMILES('[CH2]C([CH2])C(=C)CO[O]'),
    E0 = (397.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,250.578,250.623],'cm^-1')),
        HinderedRotor(inertia=(0.15211,'amu*angstrom^2'), symmetry=1, barrier=(6.77818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00404135,'amu*angstrom^2'), symmetry=1, barrier=(6.77827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00404128,'amu*angstrom^2'), symmetry=1, barrier=(6.77826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152075,'amu*angstrom^2'), symmetry=1, barrier=(6.77823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47316,'amu*angstrom^2'), symmetry=1, barrier=(65.6453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281908,0.083073,-9.17499e-05,5.98146e-08,-1.60295e-11,47928.1,34.3289], Tmin=(100,'K'), Tmax=(904.317,'K')), NASAPolynomial(coeffs=[10.8958,0.036125,-1.3876e-05,2.40511e-09,-1.58399e-13,46008.4,-15.8094], Tmin=(904.317,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])CO[O](1076)',
    structure = SMILES('[CH2]C(=C)C([CH2])CO[O]'),
    E0 = (337.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,303.017,2399.14],'cm^-1')),
        HinderedRotor(inertia=(0.175903,'amu*angstrom^2'), symmetry=1, barrier=(11.4614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13453,'amu*angstrom^2'), symmetry=1, barrier=(73.9227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00183597,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175903,'amu*angstrom^2'), symmetry=1, barrier=(11.4613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13454,'amu*angstrom^2'), symmetry=1, barrier=(73.9227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169401,0.0831285,-8.20542e-05,4.49883e-08,-1.00765e-11,40684.9,32.1791], Tmin=(100,'K'), Tmax=(1073.82,'K')), NASAPolynomial(coeffs=[13.3512,0.0340264,-1.34646e-05,2.40561e-09,-1.6266e-13,37854,-32.3542], Tmin=(1073.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)CO[O](776)',
    structure = SMILES('[CH2]C(C=C)CO[O]'),
    E0 = (220.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.109552,'amu*angstrom^2'), symmetry=1, barrier=(2.51882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110213,'amu*angstrom^2'), symmetry=1, barrier=(2.53402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110361,'amu*angstrom^2'), symmetry=1, barrier=(2.53742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91501,'amu*angstrom^2'), symmetry=1, barrier=(21.0379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.996408,0.0658248,-6.0808e-05,3.23387e-08,-7.18193e-12,26627.3,28.8667], Tmin=(100,'K'), Tmax=(1070.48,'K')), NASAPolynomial(coeffs=[10.2291,0.0313256,-1.24665e-05,2.23304e-09,-1.51098e-13,24650.6,-16.3048], Tmin=(1070.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH2])C=C(496)',
    structure = SMILES('[CH2]C([CH2])C=C'),
    E0 = (361.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,444.303],'cm^-1')),
        HinderedRotor(inertia=(0.0434666,'amu*angstrom^2'), symmetry=1, barrier=(6.11473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.043652,'amu*angstrom^2'), symmetry=1, barrier=(6.1174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0394383,'amu*angstrom^2'), symmetry=1, barrier=(69.5439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03654,0.0357961,3.13511e-06,-3.00403e-08,1.51068e-11,43603,21.9965], Tmin=(100,'K'), Tmax=(906.7,'K')), NASAPolynomial(coeffs=[9.64737,0.0219243,-6.51387e-06,1.02236e-09,-6.65462e-14,41412.9,-18.4424], Tmin=(906.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH2])C([CH2])[CH2](485)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[CH2]'),
    E0 = (617.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,711.401],'cm^-1')),
        HinderedRotor(inertia=(0.0740442,'amu*angstrom^2'), symmetry=1, barrier=(1.70242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510071,'amu*angstrom^2'), symmetry=1, barrier=(1.83404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0708634,'amu*angstrom^2'), symmetry=1, barrier=(1.62929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768939,'amu*angstrom^2'), symmetry=1, barrier=(1.76794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00487409,'amu*angstrom^2'), symmetry=1, barrier=(1.74828,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.64,'J/mol'), sigma=(6.11673,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.43 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6689,0.0650717,-5.68349e-05,2.94125e-08,-6.08673e-12,74344,30.5296], Tmin=(100,'K'), Tmax=(1324.25,'K')), NASAPolynomial(coeffs=[11.0931,0.0274373,-7.24249e-06,9.40708e-10,-4.98543e-14,72122.2,-20.6535], Tmin=(1324.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH2](499)',
    structure = SMILES('[CH2][CH]C([CH2])[CH2]'),
    E0 = (636.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1200.87],'cm^-1')),
        HinderedRotor(inertia=(0.112344,'amu*angstrom^2'), symmetry=1, barrier=(2.58301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112403,'amu*angstrom^2'), symmetry=1, barrier=(2.58438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252635,'amu*angstrom^2'), symmetry=1, barrier=(2.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252511,'amu*angstrom^2'), symmetry=1, barrier=(2.5843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95305,0.0434597,-3.02543e-05,1.30344e-08,-2.49244e-12,76589.1,26.1277], Tmin=(100,'K'), Tmax=(1189.32,'K')), NASAPolynomial(coeffs=[6.71646,0.027439,-1.00484e-05,1.70807e-09,-1.11583e-13,75456,2.32115], Tmin=(1189.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](CO[O])C([CH2])[CH2](1077)',
    structure = SMILES('[CH2][C](CO[O])C([CH2])[CH2]'),
    E0 = (666.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,874.683,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168142,'amu*angstrom^2'), symmetry=1, barrier=(7.28338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168142,'amu*angstrom^2'), symmetry=1, barrier=(7.28338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168142,'amu*angstrom^2'), symmetry=1, barrier=(7.28338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168142,'amu*angstrom^2'), symmetry=1, barrier=(7.28338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168142,'amu*angstrom^2'), symmetry=1, barrier=(7.28338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168142,'amu*angstrom^2'), symmetry=1, barrier=(7.28338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881645,0.0729937,-6.17654e-05,1.70317e-08,7.26726e-12,80280.3,38.4972], Tmin=(100,'K'), Tmax=(653.79,'K')), NASAPolynomial(coeffs=[9.04665,0.0376014,-1.39756e-05,2.36996e-09,-1.53867e-13,78901.4,0.195413], Tmin=(653.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C([CH2])CO[O](1078)',
    structure = SMILES('[CH2][C]([CH2])C([CH2])CO[O]'),
    E0 = (661.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,398.282,3093.23],'cm^-1')),
        HinderedRotor(inertia=(0.0434019,'amu*angstrom^2'), symmetry=1, barrier=(4.26497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0434019,'amu*angstrom^2'), symmetry=1, barrier=(4.26497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0434019,'amu*angstrom^2'), symmetry=1, barrier=(4.26497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0434019,'amu*angstrom^2'), symmetry=1, barrier=(4.26497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0434019,'amu*angstrom^2'), symmetry=1, barrier=(4.26497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0434019,'amu*angstrom^2'), symmetry=1, barrier=(4.26497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0186662,0.0988486,-0.000154388,1.35973e-07,-4.55255e-11,79644.6,37.8348], Tmin=(100,'K'), Tmax=(919.983,'K')), NASAPolynomial(coeffs=[6.11535,0.0436692,-1.76723e-05,3.02522e-09,-1.92231e-13,79736.1,15.5246], Tmin=(919.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])[CH]O[O](1079)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[CH]O[O]'),
    E0 = (664.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1600],'cm^-1')),
        HinderedRotor(inertia=(0.151003,'amu*angstrom^2'), symmetry=1, barrier=(3.47186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151003,'amu*angstrom^2'), symmetry=1, barrier=(3.47186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151003,'amu*angstrom^2'), symmetry=1, barrier=(3.47186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151003,'amu*angstrom^2'), symmetry=1, barrier=(3.47186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151003,'amu*angstrom^2'), symmetry=1, barrier=(3.47186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151003,'amu*angstrom^2'), symmetry=1, barrier=(3.47186,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.109194,0.0947043,-0.000126204,9.62982e-08,-2.89525e-11,80035.7,37.3078], Tmin=(100,'K'), Tmax=(936.725,'K')), NASAPolynomial(coeffs=[10.6473,0.0365739,-1.35853e-05,2.24592e-09,-1.4096e-13,78555.7,-11.0265], Tmin=(936.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[C](C)CO[O](1080)',
    structure = SMILES('[CH2]C([CH2])[C](C)CO[O]'),
    E0 = (461.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.629661,0.0748695,-6.64621e-05,3.62199e-08,-8.49227e-12,55626.9,37.5743], Tmin=(100,'K'), Tmax=(1004.38,'K')), NASAPolynomial(coeffs=[9.07521,0.0412348,-1.62301e-05,2.878e-09,-1.93183e-13,53930.4,-3.20757], Tmin=(1004.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C([CH2])CO[O](831)',
    structure = SMILES('[CH2][C](C)C([CH2])CO[O]'),
    E0 = (456.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,214.038,1757.22],'cm^-1')),
        HinderedRotor(inertia=(0.0845112,'amu*angstrom^2'), symmetry=1, barrier=(2.8343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845112,'amu*angstrom^2'), symmetry=1, barrier=(2.8343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845112,'amu*angstrom^2'), symmetry=1, barrier=(2.8343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845112,'amu*angstrom^2'), symmetry=1, barrier=(2.8343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845112,'amu*angstrom^2'), symmetry=1, barrier=(2.8343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0845112,'amu*angstrom^2'), symmetry=1, barrier=(2.8343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.041255,0.0966887,-0.00014029,1.21599e-07,-4.12528e-11,54979.5,36.6654], Tmin=(100,'K'), Tmax=(889.565,'K')), NASAPolynomial(coeffs=[6.07295,0.0474729,-2.0047e-05,3.56561e-09,-2.34485e-13,54780.5,13.1853], Tmin=(889.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])[CH]OO(1081)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[CH]OO'),
    E0 = (512.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0123932,0.0929842,-9.85622e-05,4.90145e-08,-4.0724e-12,61748.9,36.1577], Tmin=(100,'K'), Tmax=(677.84,'K')), NASAPolynomial(coeffs=[12.0081,0.038677,-1.48557e-05,2.55692e-09,-1.67189e-13,59744,-19.843], Tmin=(677.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C(C)CO[O](1082)',
    structure = SMILES('[CH2][C]([CH2])C(C)CO[O]'),
    E0 = (456.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.041255,0.0966887,-0.00014029,1.21599e-07,-4.12528e-11,54979.5,35.9723], Tmin=(100,'K'), Tmax=(889.565,'K')), NASAPolynomial(coeffs=[6.07295,0.0474729,-2.0047e-05,3.56561e-09,-2.34485e-13,54780.5,12.4921], Tmin=(889.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C(C)[CH]O[O](1083)',
    structure = SMILES('[CH2]C([CH2])C(C)[CH]O[O]'),
    E0 = (459.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.038573,0.0919388,-0.000109796,7.86364e-08,-2.31377e-11,55368.6,35.2755], Tmin=(100,'K'), Tmax=(865.468,'K')), NASAPolynomial(coeffs=[10.4805,0.0405974,-1.60901e-05,2.81764e-09,-1.85851e-13,53649.9,-13.3634], Tmin=(865.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](CO[O])C([CH2])C(830)',
    structure = SMILES('[CH2][C](CO[O])C([CH2])C'),
    E0 = (461.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,299.276,2980.39],'cm^-1')),
        HinderedRotor(inertia=(0.0281756,'amu*angstrom^2'), symmetry=1, barrier=(3.71972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281756,'amu*angstrom^2'), symmetry=1, barrier=(3.71972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281756,'amu*angstrom^2'), symmetry=1, barrier=(3.71972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281756,'amu*angstrom^2'), symmetry=1, barrier=(3.71972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281756,'amu*angstrom^2'), symmetry=1, barrier=(3.71972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281756,'amu*angstrom^2'), symmetry=1, barrier=(3.71972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.629661,0.0748695,-6.64621e-05,3.62199e-08,-8.49227e-12,55626.9,38.2674], Tmin=(100,'K'), Tmax=(1004.38,'K')), NASAPolynomial(coeffs=[9.07521,0.0412348,-1.62301e-05,2.878e-09,-1.93183e-13,53930.4,-2.51442], Tmin=(1004.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C([CH2])[CH]O[O](832)',
    structure = SMILES('[CH2]C(C)C([CH2])[CH]O[O]'),
    E0 = (459.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.038573,0.0919388,-0.000109796,7.86364e-08,-2.31377e-11,55368.6,35.9687], Tmin=(100,'K'), Tmax=(865.468,'K')), NASAPolynomial(coeffs=[10.4805,0.0405974,-1.60901e-05,2.81764e-09,-1.85851e-13,53649.9,-12.6702], Tmin=(865.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](COO)C([CH2])[CH2](1084)',
    structure = SMILES('[CH2][C](COO)C([CH2])[CH2]'),
    E0 = (514.582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,218.508,1955.54],'cm^-1')),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128136,'amu*angstrom^2'), symmetry=1, barrier=(3.59915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277007,0.0814272,-7.88134e-05,4.53725e-08,-1.08896e-11,62024.4,39.8607], Tmin=(100,'K'), Tmax=(999.063,'K')), NASAPolynomial(coeffs=[10.9614,0.0386505,-1.45895e-05,2.51716e-09,-1.65953e-13,59889.5,-11.6753], Tmin=(999.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C([CH2])COO(1085)',
    structure = SMILES('[CH2][C]([CH2])C([CH2])COO'),
    E0 = (509.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1600],'cm^-1')),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151358,'amu*angstrom^2'), symmetry=1, barrier=(3.48002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.240818,0.102377,-0.000149413,1.26258e-07,-4.15786e-11,61373.9,38.0079], Tmin=(100,'K'), Tmax=(903.002,'K')), NASAPolynomial(coeffs=[7.82439,0.0451217,-1.85422e-05,3.23712e-09,-2.09954e-13,60795.1,4.78099], Tmin=(903.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(403.252,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (475.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (635.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (635.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (632.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (720.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (938.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (932.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (930.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (930.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (483.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (483.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (482.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (484.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (558.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (537.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (615.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (498.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (539.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (539.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (500.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (500.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (609.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (560.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (602.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (570.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (603.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (534.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (608.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (841.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (842.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (878.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (872.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (876.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (578.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (617.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (641.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (593.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (592.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (593.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (538.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (598.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (550.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['allyl(82)', 'allylperoxy(462)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])C[CH]CO[O](743)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])[CH]CCO[O](1061)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2][CH]CC([CH2])CO[O](742)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2]C([CH2])C([CH2])C[O](1062)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2]C([CH2])[CH]CO[O](1063)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])CO[O](778)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(CO[O])C([CH2])[CH2](1064)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([CH2])C([CH2])CO[O](1065)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C(CO[O])C1CC1(1066)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C1CCC1CO[O](1067)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])C1COOC1(1068)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C1COOCC1[CH2](1069)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.42e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R6_SSSSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSSSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['O(4)', '[CH2]C([CH2])C1COC1(1070)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['O(4)', '[CH2]C1COCC1[CH2](1071)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.026e+11,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO_SSS;C_pri_rad_intra;OO] for rate rule [R4OO_SSS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['HO2(9)', '[CH2]C(=C)C([CH2])[CH2](635)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C(CO[O])C(=C)C(829)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C(C)C(=C)CO[O](828)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C(=C)C(C)CO[O](1072)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])C(=C)COO(1073)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C(=C)C([CH2])COO(1074)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C([CH2])C(=C)CO[O](1075)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C(=C)C([CH2])CO[O](1076)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(20)', '[CH2]C(C=C)CO[O](776)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH2](497)', 'allylperoxy(462)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.540096,'m^3/(mol*s)'), n=2.05449, Ea=(13.6169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;CJ] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]O[O](46)', '[CH2]C([CH2])C=C(496)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00140605,'m^3/(mol*s)'), n=2.41796, Ea=(36.1629,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CsJ] for rate rule [Cds-CsH_Cds-HH;CsJ-OsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['allyl(82)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O2(2)', '[CH2]C([CH2])C([CH2])[CH2](485)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.1851e+07,'m^3/(mol*s)'), n=-0.119415, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 8.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(308.407,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]O[O](46)', '[CH2][CH]C([CH2])[CH2](499)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][C](CO[O])C([CH2])[CH2](1077)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][C]([CH2])C([CH2])CO[O](1078)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2]C([CH2])C([CH2])[CH]O[O](1079)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])[C](C)CO[O](1080)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2][C](C)C([CH2])CO[O](831)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])C([CH2])[CH]OO(1081)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 242 used for R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2][C]([CH2])C(C)CO[O](1082)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C([CH2])C(C)[CH]O[O](1083)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2][C](CO[O])C([CH2])C(830)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    products = ['[CH2]C(C)C([CH2])[CH]O[O](832)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(608,'s^-1'), n=2.77, Ea=(62.3834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 324 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][C](COO)C([CH2])[CH2](1084)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][C]([CH2])C([CH2])COO(1085)'],
    products = ['[CH2]C([CH2])C([CH2])CO[O](744)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #314',
    isomers = [
        '[CH2]C([CH2])C([CH2])CO[O](744)',
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
    label = 'PDepNetwork #314',
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


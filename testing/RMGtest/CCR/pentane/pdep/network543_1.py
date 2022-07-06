species(
    label = '[CH2][CH]OC([CH2])CO[O](1494)',
    structure = SMILES('[CH2][CH]OC([CH2])CO[O]'),
    E0 = (358.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.981562,0.119585,-0.000196359,1.66279e-07,-5.41663e-11,43291.5,34.8217], Tmin=(100,'K'), Tmax=(877.105,'K')), NASAPolynomial(coeffs=[13.4791,0.0337002,-1.53861e-05,2.80994e-09,-1.8649e-13,41521.6,-28.6749], Tmin=(877.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJC(C)OC) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = 'vinoxy(1351)',
    structure = SMILES('[CH2]C=O'),
    E0 = (6.53237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1127.76,1127.76,1127.77,1945],'cm^-1')),
        HinderedRotor(inertia=(0.00342653,'amu*angstrom^2'), symmetry=1, barrier=(3.09254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51396,0.0045106,2.69938e-05,-3.71579e-08,1.42465e-11,808.756,8.86322], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[6.49564,0.00770518,-2.52913e-06,4.69083e-10,-3.5128e-14,-479.661,-9.13857], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.53237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(CO[O])OC=C(1506)',
    structure = SMILES('[CH2]C(CO[O])OC=C'),
    E0 = (52.2622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.772442,0.0982834,-0.000112638,6.43146e-08,-1.42085e-11,6463.3,31.6881], Tmin=(100,'K'), Tmax=(1117.65,'K')), NASAPolynomial(coeffs=[21.2811,0.0193544,-6.70616e-06,1.12688e-09,-7.42904e-14,1533.71,-77.1602], Tmin=(1117.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH]O[CH]CCO[O](2106)',
    structure = SMILES('[CH2][CH]O[CH]CCO[O]'),
    E0 = (339.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.912311,0.115818,-0.000182232,1.48891e-07,-4.72051e-11,40945.4,34.612], Tmin=(100,'K'), Tmax=(872.578,'K')), NASAPolynomial(coeffs=[14.516,0.0315238,-1.40015e-05,2.53872e-09,-1.68406e-13,38769.5,-34.758], Tmin=(872.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOCs) + radical(CCsJOCs) + radical(CJCO)"""),
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
    label = '[CH2]C([O])CO[O](1369)',
    structure = SMILES('[CH2]C([O])CO[O]'),
    E0 = (212.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,203.269,203.415],'cm^-1')),
        HinderedRotor(inertia=(0.392719,'amu*angstrom^2'), symmetry=1, barrier=(11.5392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00489682,'amu*angstrom^2'), symmetry=1, barrier=(48.6712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393606,'amu*angstrom^2'), symmetry=1, barrier=(11.5369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (89.07,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28256,0.0633366,-8.80423e-05,6.82761e-08,-2.13442e-11,25668.9,24.1399], Tmin=(100,'K'), Tmax=(801.043,'K')), NASAPolynomial(coeffs=[9.19638,0.022794,-1.02044e-05,1.89819e-09,-1.29584e-13,24433.9,-12.0788], Tmin=(801.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO)"""),
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
    label = '[CH2][CH]OC([CH2])C[O](2107)',
    structure = SMILES('[CH2][CH]OC([CH2])C[O]'),
    E0 = (360.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196959,0.103381,-0.000172394,1.532e-07,-5.2337e-11,43526.2,31.345], Tmin=(100,'K'), Tmax=(867.036,'K')), NASAPolynomial(coeffs=[9.18224,0.0369572,-1.74239e-05,3.24336e-09,-2.18212e-13,42770,-7.54739], Tmin=(867.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOCs) + radical(CJCO)"""),
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
    label = '[CH2][CH]O[CH]CO[O](2108)',
    structure = SMILES('[CH2][CH]O[CH]CO[O]'),
    E0 = (362.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.2302,0.100601,-0.000166767,1.38487e-07,-4.39197e-11,43781.2,29.9188], Tmin=(100,'K'), Tmax=(893.998,'K')), NASAPolynomial(coeffs=[13.4736,0.0233193,-1.03077e-05,1.83495e-09,-1.19145e-13,41969,-31.09], Tmin=(893.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOCs) + radical(CCsJOCs) + radical(CJCO)"""),
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
    label = '[CH2][C]OC([CH2])CO[O](2109)',
    structure = SMILES('[CH2][C]OC([CH2])CO[O]'),
    E0 = (639.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.870658,0.118254,-0.000201475,1.74798e-07,-5.81576e-11,77026.6,33.7096], Tmin=(100,'K'), Tmax=(866.065,'K')), NASAPolynomial(coeffs=[12.9289,0.0328036,-1.58652e-05,2.96917e-09,-1.99784e-13,75450.8,-26.1789], Tmin=(866.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJC(C)OC) + radical(CJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C(CO[O])O[CH][CH2](2110)',
    structure = SMILES('[CH]C(CO[O])O[CH][CH2]'),
    E0 = (596.238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.6508,0.109496,-0.000171012,1.39425e-07,-4.44013e-11,71871.5,34.8871], Tmin=(100,'K'), Tmax=(853.82,'K')), NASAPolynomial(coeffs=[14.0179,0.0304261,-1.39187e-05,2.56903e-09,-1.72748e-13,69743.9,-31.3535], Tmin=(853.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOCs) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]OC([CH2])CO[O](2111)',
    structure = SMILES('[CH][CH]OC([CH2])CO[O]'),
    E0 = (595.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.88166,0.11745,-0.000196119,1.66777e-07,-5.4371e-11,71747.3,34.8334], Tmin=(100,'K'), Tmax=(878.453,'K')), NASAPolynomial(coeffs=[13.553,0.0314168,-1.45397e-05,2.66132e-09,-1.76497e-13,69994.8,-28.4755], Tmin=(878.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJC(C)OC) + radical(CCsJOCs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC(CO[O])O1(1495)',
    structure = SMILES('[CH2]C1CC(CO[O])O1'),
    E0 = (105.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465936,0.074749,-6.11956e-05,1.73009e-08,3.32186e-12,12762.1,27.251], Tmin=(100,'K'), Tmax=(817.395,'K')), NASAPolynomial(coeffs=[13.8185,0.0275945,-8.03835e-06,1.16728e-09,-6.93183e-14,9971.68,-38.1918], Tmin=(817.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH]OC1COOC1(2112)',
    structure = SMILES('[CH2][CH]OC1COOC1'),
    E0 = (78.3524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.746801,0.0850326,-8.22509e-05,4.10149e-08,-7.7444e-12,9611.12,30.9015], Tmin=(100,'K'), Tmax=(1502.72,'K')), NASAPolynomial(coeffs=[19.2769,0.0187878,-3.20425e-06,2.14047e-10,-2.85979e-15,5054.73,-68.993], Tmin=(1502.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.3524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1COOC([CH2])O1(2113)',
    structure = SMILES('[CH2]C1COOC([CH2])O1'),
    E0 = (60.5879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68613,0.0886105,-8.4689e-05,4.0928e-08,-7.22805e-12,7523.37,33.9035], Tmin=(100,'K'), Tmax=(1717.33,'K')), NASAPolynomial(coeffs=[18.3983,0.0166678,1.27064e-07,-5.29338e-10,4.91098e-14,4335.48,-63.0504], Tmin=(1717.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.5879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(124trioxane) + radical(CJC(C)OC) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][CH]OC1COC1(2114)',
    structure = SMILES('[CH2][CH]OC1COC1'),
    E0 = (106.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0659716,0.072702,-6.84458e-05,3.47131e-08,-6.76607e-12,12956.9,26.4998], Tmin=(100,'K'), Tmax=(1442.69,'K')), NASAPolynomial(coeffs=[15.2993,0.0207585,-4.34581e-06,4.28464e-10,-1.67536e-14,9571.71,-49.0741], Tmin=(1442.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1COC([CH2])O1(2115)',
    structure = SMILES('[CH2]C1COC([CH2])O1'),
    E0 = (13.8545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.416044,0.0736332,-7.73402e-05,4.11015e-08,-7.91071e-12,1845.57,26.058], Tmin=(100,'K'), Tmax=(1580.57,'K')), NASAPolynomial(coeffs=[15.749,0.0114338,1.89335e-06,-8.40453e-10,7.12195e-14,-605.073,-50.9165], Tmin=(1580.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.8545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(1,3-Dioxolane) + radical(CJC(C)OC) + radical(CJCO)"""),
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
    label = '[CH2][CH]OC([CH2])=C(1863)',
    structure = SMILES('[CH2][CH]OC([CH2])=C'),
    E0 = (337.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,409.215,409.221,409.253],'cm^-1')),
        HinderedRotor(inertia=(0.135313,'amu*angstrom^2'), symmetry=1, barrier=(16.0803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135316,'amu*angstrom^2'), symmetry=1, barrier=(16.0806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135329,'amu*angstrom^2'), symmetry=1, barrier=(16.0804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135307,'amu*angstrom^2'), symmetry=1, barrier=(16.0804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0881852,0.0768561,-8.49605e-05,4.51228e-08,-9.02004e-12,40706.3,26.1587], Tmin=(100,'K'), Tmax=(1359.75,'K')), NASAPolynomial(coeffs=[20.7136,0.00828131,-1.16948e-06,4.87588e-11,1.20998e-15,35731.7,-78.0811], Tmin=(1359.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=C(O)CJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]COC(=C)CO[O](2116)',
    structure = SMILES('[CH2]COC(=C)CO[O]'),
    E0 = (54.6645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00212011,0.0871228,-9.4373e-05,5.36242e-08,-1.21678e-11,6719.35,30.9512], Tmin=(100,'K'), Tmax=(1070.25,'K')), NASAPolynomial(coeffs=[15.7687,0.0281958,-1.17843e-05,2.17898e-09,-1.50689e-13,3344.51,-46.184], Tmin=(1070.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.6645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]OC(=C)COO(2117)',
    structure = SMILES('[CH2][CH]OC(=C)COO'),
    E0 = (96.5872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.681709,0.100731,-0.0001187,6.97023e-08,-1.59776e-11,11787.4,32.8132], Tmin=(100,'K'), Tmax=(1070.99,'K')), NASAPolynomial(coeffs=[20.175,0.0228343,-9.59977e-06,1.78953e-09,-1.24761e-13,7319.93,-69.2386], Tmin=(1070.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.5872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]OC(=C)CO[O](2118)',
    structure = SMILES('[CH2][CH]OC(=C)CO[O]'),
    E0 = (248.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.407234,0.097001,-0.000122888,7.83585e-08,-1.9493e-11,30057.4,32.5883], Tmin=(100,'K'), Tmax=(989.586,'K')), NASAPolynomial(coeffs=[18.1756,0.021885,-9.02429e-06,1.64797e-09,-1.12923e-13,26379.7,-56.8678], Tmin=(989.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCsJOC(O)) + radical(CJCO)"""),
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
    label = '[CH2][CH]OC=C(1522)',
    structure = SMILES('[CH2][CH]OC=C'),
    E0 = (220.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.952254,'amu*angstrom^2'), symmetry=1, barrier=(21.8942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952621,'amu*angstrom^2'), symmetry=1, barrier=(21.9026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.952086,'amu*angstrom^2'), symmetry=1, barrier=(21.8903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03469,0.0521468,-2.34313e-05,-2.40782e-08,1.80817e-11,26580,19.4356], Tmin=(100,'K'), Tmax=(911.623,'K')), NASAPolynomial(coeffs=[20.5653,0.0015492,2.07159e-06,-4.94943e-10,3.23855e-14,21560.7,-80.9802], Tmin=(911.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH][O](1556)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (367.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1667.93],'cm^-1')),
        HinderedRotor(inertia=(0.00517725,'amu*angstrom^2'), symmetry=1, barrier=(10.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22654,0.0212781,-3.59354e-05,3.91027e-08,-1.59281e-11,44278.5,12.2199], Tmin=(100,'K'), Tmax=(830.699,'K')), NASAPolynomial(coeffs=[2.17156,0.016018,-7.76586e-06,1.51127e-09,-1.05387e-13,44810.5,19.2612], Tmin=(830.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
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
    label = '[CH2][CH]OC([CH2])[CH2](1479)',
    structure = SMILES('[CH2][CH]OC([CH2])[CH2]'),
    E0 = (505.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3228.63,'J/mol'), sigma=(5.97178,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.30 K, Pc=34.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.208136,0.0993995,-0.000158986,1.30207e-07,-4.10134e-11,60926,26.8809], Tmin=(100,'K'), Tmax=(891.629,'K')), NASAPolynomial(coeffs=[13.1919,0.0253305,-1.09034e-05,1.93371e-09,-1.25921e-13,59091.1,-33.1188], Tmin=(891.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]O[CH][CH2](1569)',
    structure = SMILES('[CH2][CH]O[CH][CH2]'),
    E0 = (510.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,278.26,278.282,278.332],'cm^-1')),
        HinderedRotor(inertia=(0.00217771,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0021765,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296218,'amu*angstrom^2'), symmetry=1, barrier=(16.2788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296261,'amu*angstrom^2'), symmetry=1, barrier=(16.278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866249,0.0712249,-9.91396e-05,6.6928e-08,-1.6514e-11,61536,21.71], Tmin=(100,'K'), Tmax=(766.355,'K')), NASAPolynomial(coeffs=[13.5583,0.0141332,-5.31156e-06,8.93057e-10,-5.71145e-14,59321.8,-37.8983], Tmin=(766.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOCs) + radical(CCsJOCs) + radical(CJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]O[C]([CH2])CO[O](2119)',
    structure = SMILES('[CH2][CH]O[C]([CH2])CO[O]'),
    E0 = (539.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09192,0.130444,-0.000247655,2.28642e-07,-7.82098e-11,65021.2,35.163], Tmin=(100,'K'), Tmax=(897.743,'K')), NASAPolynomial(coeffs=[10.0051,0.0369266,-1.77604e-05,3.23609e-09,-2.10885e-13,64804.8,-7.28519], Tmin=(897.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CCsJOCs) + radical(CJC(C)OC) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]OC([CH2])[CH]O[O](2120)',
    structure = SMILES('[CH2][CH]OC([CH2])[CH]O[O]'),
    E0 = (547.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.998606,0.122572,-0.00021407,1.87407e-07,-6.21528e-11,65970.9,35.0318], Tmin=(100,'K'), Tmax=(887.302,'K')), NASAPolynomial(coeffs=[12.7185,0.0329051,-1.54408e-05,2.8204e-09,-1.85779e-13,64632.1,-23.3319], Tmin=(887.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCsJOCs) + radical(CJC(C)OC) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]O[C](C)CO[O](2121)',
    structure = SMILES('[CH2][CH]O[C](C)CO[O]'),
    E0 = (328.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.712172,0.117353,-0.000203043,1.82554e-07,-6.23204e-11,39693.4,34.0477], Tmin=(100,'K'), Tmax=(879.854,'K')), NASAPolynomial(coeffs=[9.6522,0.0397215,-1.86759e-05,3.44415e-09,-2.29302e-13,39050.6,-7.91607], Tmin=(879.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = '[CH2]CO[C]([CH2])CO[O](2122)',
    structure = SMILES('[CH2]CO[C]([CH2])CO[O]'),
    E0 = (358.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.773086,0.123231,-0.000227892,2.14374e-07,-7.52279e-11,43305.9,34.5724], Tmin=(100,'K'), Tmax=(884.311,'K')), NASAPolynomial(coeffs=[7.1943,0.0443624,-2.1462e-05,3.98042e-09,-2.64705e-13,43571.5,6.58269], Tmin=(884.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH]OC([CH2])[CH]OO(2123)',
    structure = SMILES('[CH2][CH]OC([CH2])[CH]OO'),
    E0 = (395.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24954,0.126006,-0.00020881,1.77433e-07,-5.81663e-11,47699.9,35.174], Tmin=(100,'K'), Tmax=(863.646,'K')), NASAPolynomial(coeffs=[14.3441,0.0345015,-1.63945e-05,3.05224e-09,-2.05165e-13,45725.5,-33.607], Tmin=(863.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(CJC(C)OC) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]OC(C)[CH]O[O](2124)',
    structure = SMILES('[CH2][CH]OC(C)[CH]O[O]'),
    E0 = (336.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.596603,0.109207,-0.000168459,1.39988e-07,-4.56934e-11,40642.1,33.8376], Tmin=(100,'K'), Tmax=(849.911,'K')), NASAPolynomial(coeffs=[12.2685,0.03587,-1.64565e-05,3.05247e-09,-2.0621e-13,38917.2,-23.4196], Tmin=(849.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](CO[O])O[CH]C(2125)',
    structure = SMILES('[CH2][C](CO[O])O[CH]C'),
    E0 = (327.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.927302,0.125111,-0.000227418,2.08898e-07,-7.1839e-11,39568.5,33.9383], Tmin=(100,'K'), Tmax=(889.387,'K')), NASAPolynomial(coeffs=[9.12911,0.0408146,-1.93573e-05,3.55097e-09,-2.34272e-13,39324.9,-4.71261], Tmin=(889.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CCsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]COC([CH2])[CH]O[O](2126)',
    structure = SMILES('[CH2]COC([CH2])[CH]O[O]'),
    E0 = (366.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.672876,0.115278,-0.00019403,1.72806e-07,-5.90511e-11,44255.3,34.4166], Tmin=(100,'K'), Tmax=(867.786,'K')), NASAPolynomial(coeffs=[9.8652,0.0404147,-1.91856e-05,3.57502e-09,-2.40461e-13,43416.2,-9.226], Tmin=(867.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)OC) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]O[C]([CH2])COO(2127)',
    structure = SMILES('[CH2][CH]O[C]([CH2])COO'),
    E0 = (387.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35522,0.134028,-0.000242936,2.19375e-07,-7.45174e-11,46750.7,35.3491], Tmin=(100,'K'), Tmax=(884.162,'K')), NASAPolynomial(coeffs=[11.6892,0.0384208,-1.86541e-05,3.45357e-09,-2.29067e-13,45874.4,-17.8882], Tmin=(884.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOCs) + radical(CCsJOCs) + radical(CJC(C)OC) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH]O[O])O[CH]C(2128)',
    structure = SMILES('[CH2]C([CH]O[O])O[CH]C'),
    E0 = (335.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.827336,0.117158,-0.000193544,1.6729e-07,-5.56283e-11,40517.9,33.7835], Tmin=(100,'K'), Tmax=(873.116,'K')), NASAPolynomial(coeffs=[11.81,0.0368496,-1.70709e-05,3.1432e-09,-2.0983e-13,39165.4,-20.5774], Tmin=(873.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCsJOCs) + radical(CJC(C)OC)"""),
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
    E0 = (358.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (358.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (517.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (770.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (603.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (800.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (850.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (808.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (806.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (366.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (365.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (367.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (441.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (451.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (497.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (421.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (383.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (471.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (394.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (458.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (724.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (496.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (716.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (751.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (758.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (461.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (523.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (524.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (475.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (489.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (463.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (471.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (413.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['vinoxy(1351)', 'allylperoxy(462)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2]C(CO[O])OC=C(1506)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][CH]O[CH]CCO[O](2106)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][CH2](502)', '[CH2]C([O])CO[O](1369)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2][CH]OC([CH2])C[O](2107)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH]O[CH]CO[O](2108)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C]OC([CH2])CO[O](2109)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(CO[O])O[CH][CH2](2110)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH][CH]OC([CH2])CO[O](2111)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2]C1CC(CO[O])O1(1495)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][CH]OC1COOC1(2112)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2]C1COOC([CH2])O1(2113)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['O(4)', '[CH2][CH]OC1COC1(2114)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['O(4)', '[CH2]C1COC([CH2])O1(2115)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.14131e+09,'s^-1'), n=0.54, Ea=(93.1987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R4OO;Y_rad_intra;OO] for rate rule [R4OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['HO2(9)', '[CH2][CH]OC([CH2])=C(1863)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2]COC(=C)CO[O](2116)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][CH]OC(=C)COO(2117)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][CH]OC(=C)CO[O](2118)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['vinoxy(1351)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]O[O](46)', '[CH2][CH]OC=C(1522)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(364.868,'m^3/(mol*s)'), n=0.88, Ea=(33.1163,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CsJ-OsHH] for rate rule [Cds-OsH_Cds;CsJ-OsHH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(308.407,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O2(2)', '[CH2][CH]OC([CH2])[CH2](1479)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1233.63,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]O[O](46)', '[CH2][CH]O[CH][CH2](1569)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(203262,'m^3/(mol*s)'), n=0.35323, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0148472730165, var=2.75207767881, Tref=1000.0, N=9, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R
    Total Standard Deviation in ln(k): 3.36303735057
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH]O[C]([CH2])CO[O](2119)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH]OC([CH2])[CH]O[O](2120)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][CH]O[C](C)CO[O](2121)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2]CO[C]([CH2])CO[O](2122)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.9172e+08,'s^-1'), n=1.32036, Ea=(164.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][CH]OC([CH2])[CH]OO(2123)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.66e+08,'s^-1'), n=1.28, Ea=(166.272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 241 used for R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][CH]OC(C)[CH]O[O](2124)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2][C](CO[O])O[CH]C(2125)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]COC([CH2])[CH]O[O](2126)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(55800,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]O[C]([CH2])COO(2127)'],
    products = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]OC([CH2])CO[O](1494)'],
    products = ['[CH2]C([CH]O[O])O[CH]C(2128)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(126713,'s^-1'), n=1.75034, Ea=(55.0506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_1H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #543',
    isomers = [
        '[CH2][CH]OC([CH2])CO[O](1494)',
    ],
    reactants = [
        ('vinoxy(1351)', 'allylperoxy(462)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #543',
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


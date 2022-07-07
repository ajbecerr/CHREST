species(
    label = '[CH2]C([O])C[C]=CO[O](21021)',
    structure = SMILES('[CH2]C([O])C[C]=CO[O]'),
    E0 = (534.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,339.822,339.822,339.822],'cm^-1')),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28743,0.0846749,-0.00010982,7.56961e-08,-2.07614e-11,64443.7,33.5907], Tmin=(100,'K'), Tmax=(892.541,'K')), NASAPolynomial(coeffs=[13.3715,0.0260365,-1.12712e-05,2.08583e-09,-1.42994e-13,62108.1,-28.045], Tmin=(892.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cds_S)"""),
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
    label = '[O][CH]CC[C]=CO[O](21034)',
    structure = SMILES('[O][CH]CC[C]=CO[O]'),
    E0 = (511.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,249.331,249.331,249.331,2676.45],'cm^-1')),
        HinderedRotor(inertia=(0.179106,'amu*angstrom^2'), symmetry=1, barrier=(7.90113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179102,'amu*angstrom^2'), symmetry=1, barrier=(7.90113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17911,'amu*angstrom^2'), symmetry=1, barrier=(7.90114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179106,'amu*angstrom^2'), symmetry=1, barrier=(7.90115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4412.68,'J/mol'), sigma=(7.2056,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=689.25 K, Pc=26.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271548,0.092963,-0.000155628,1.4138e-07,-4.93389e-11,61633.6,33.9522], Tmin=(100,'K'), Tmax=(863.087,'K')), NASAPolynomial(coeffs=[7.28241,0.037016,-1.76313e-05,3.30244e-09,-2.2322e-13,61297,6.22209], Tmin=(863.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C(=C)[CH]O[O](20472)',
    structure = SMILES('[CH2]C([O])C(=C)[CH]O[O]'),
    E0 = (401.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,916.311],'cm^-1')),
        HinderedRotor(inertia=(0.0340544,'amu*angstrom^2'), symmetry=1, barrier=(3.07764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.929837,'amu*angstrom^2'), symmetry=1, barrier=(21.3788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359981,'amu*angstrom^2'), symmetry=1, barrier=(21.3844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359349,'amu*angstrom^2'), symmetry=1, barrier=(21.3742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0859817,0.0760865,-7.36512e-05,3.52531e-08,-6.5779e-12,48406.8,34.4714], Tmin=(100,'K'), Tmax=(1310.73,'K')), NASAPolynomial(coeffs=[19.2616,0.0175676,-6.68204e-06,1.19097e-09,-8.1103e-14,43380,-63.2283], Tmin=(1310.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(CJCO)"""),
)

species(
    label = 'O2(S)(666)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,10302.3,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CCC([CH2])[O](15960)',
    structure = SMILES('C#CCC([CH2])[O]'),
    E0 = (395.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,375.646,376.268],'cm^-1')),
        HinderedRotor(inertia=(0.0681379,'amu*angstrom^2'), symmetry=1, barrier=(6.68049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30015,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18883,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1544,0.0586267,-5.99396e-05,3.28582e-08,-7.1473e-12,47647.1,23.1773], Tmin=(100,'K'), Tmax=(1123.73,'K')), NASAPolynomial(coeffs=[12.4442,0.0184382,-6.29231e-06,1.03001e-09,-6.60986e-14,45109.9,-32.6059], Tmin=(1123.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = '[CH2]C([O])C[C]=C[O](22838)',
    structure = SMILES('[CH2]C([O])C[C]=C[O]'),
    E0 = (399.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,430.592,430.592,430.592,430.592],'cm^-1')),
        HinderedRotor(inertia=(0.123571,'amu*angstrom^2'), symmetry=1, barrier=(16.2583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123571,'amu*angstrom^2'), symmetry=1, barrier=(16.2583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123571,'amu*angstrom^2'), symmetry=1, barrier=(16.2583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519954,0.0685225,-6.19685e-05,2.30128e-08,-1.57337e-12,48203.9,29.3556], Tmin=(100,'K'), Tmax=(988.478,'K')), NASAPolynomial(coeffs=[17.519,0.0154744,-5.35529e-06,9.40589e-10,-6.54113e-14,44074.3,-56.348], Tmin=(988.478,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C[C]=CO[O](21352)',
    structure = SMILES('[CH2][CH]C[C]=CO[O]'),
    E0 = (664.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180,2434.84],'cm^-1')),
        HinderedRotor(inertia=(0.197963,'amu*angstrom^2'), symmetry=1, barrier=(4.55155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197958,'amu*angstrom^2'), symmetry=1, barrier=(4.55144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19796,'amu*angstrom^2'), symmetry=1, barrier=(4.55149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197957,'amu*angstrom^2'), symmetry=1, barrier=(4.55143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15987,0.0679238,-9.57544e-05,8.02508e-08,-2.69639e-11,80070.8,32.1736], Tmin=(100,'K'), Tmax=(850.899,'K')), NASAPolynomial(coeffs=[6.91059,0.0311414,-1.37273e-05,2.51921e-09,-1.69847e-13,79445,7.432], Tmin=(850.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S)"""),
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
    label = '[O][CH]C[C]=CO[O](22760)',
    structure = SMILES('[O][CH]C[C]=CO[O]'),
    E0 = (535.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,180,2690.78],'cm^-1')),
        HinderedRotor(inertia=(0.384927,'amu*angstrom^2'), symmetry=1, barrier=(8.85022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384944,'amu*angstrom^2'), symmetry=1, barrier=(8.85061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384869,'amu*angstrom^2'), symmetry=1, barrier=(8.8489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.960614,0.0776544,-0.000139785,1.30389e-07,-4.57495e-11,64469.1,29.2346], Tmin=(100,'K'), Tmax=(877.507,'K')), NASAPolynomial(coeffs=[6.23866,0.0288148,-1.39397e-05,2.59928e-09,-1.74014e-13,64496.9,9.89701], Tmin=(877.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[C]=CO[O](21346)',
    structure = SMILES('[C]=CO[O]'),
    E0 = (656.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.353101,'amu*angstrom^2'), symmetry=1, barrier=(8.11848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89582,0.0288275,-5.64855e-05,5.40338e-08,-1.90645e-11,78935.3,14.7462], Tmin=(100,'K'), Tmax=(877.636,'K')), NASAPolynomial(coeffs=[5.05128,0.00860812,-4.16047e-06,7.90553e-10,-5.33044e-14,78957.3,6.90957], Tmin=(877.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])[O](1366)',
    structure = SMILES('[CH2]C([CH2])[O]'),
    E0 = (362.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,450.564],'cm^-1')),
        HinderedRotor(inertia=(0.000831827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.050104,'amu*angstrom^2'), symmetry=1, barrier=(7.38455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17517,0.0359905,-3.15384e-05,1.43889e-08,-2.59785e-12,43692.5,16.9615], Tmin=(100,'K'), Tmax=(1341.41,'K')), NASAPolynomial(coeffs=[10.4304,0.0113739,-4.01123e-06,7.08062e-10,-4.81153e-14,41477.8,-25.2898], Tmin=(1341.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CC(C)OJ) + radical(CJCO) + radical(CJCO)"""),
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
    label = '[CH]C([O])C[C]=CO[O](21141)',
    structure = SMILES('[CH]C([O])C[C]=CO[O]'),
    E0 = (771.352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381851,0.0826087,-0.000109836,7.65572e-08,-2.11363e-11,92899.8,33.6218], Tmin=(100,'K'), Tmax=(888.245,'K')), NASAPolynomial(coeffs=[13.4543,0.0237377,-1.04158e-05,1.93505e-09,-1.32821e-13,90577.6,-27.8957], Tmin=(888.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]OC=[C]CC1CO1(22839)',
    structure = SMILES('[O]OC=[C]CC1CO1'),
    E0 = (282.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.260359,0.0746919,-8.17045e-05,4.77791e-08,-1.06905e-11,34098.8,30.2673], Tmin=(100,'K'), Tmax=(1250.03,'K')), NASAPolynomial(coeffs=[14.6459,0.0193083,-4.02452e-06,3.66302e-10,-1.13073e-14,31233,-39.4224], Tmin=(1250.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])CC1=COO1(22840)',
    structure = SMILES('[CH2]C([O])CC1=COO1'),
    E0 = (325.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.555542,0.0635701,-3.59994e-05,-4.57656e-09,7.51004e-12,39244.3,30.6132], Tmin=(100,'K'), Tmax=(1012.23,'K')), NASAPolynomial(coeffs=[18.3202,0.0188167,-7.38869e-06,1.41502e-09,-1.0343e-13,34344.2,-61.7462], Tmin=(1012.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1CC(=CO[O])O1(22821)',
    structure = SMILES('[CH2]C1CC(=CO[O])O1'),
    E0 = (186.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00173,0.0828597,-8.61923e-05,4.30276e-08,-7.85451e-12,22641,29.8843], Tmin=(100,'K'), Tmax=(1615.69,'K')), NASAPolynomial(coeffs=[21.4156,0.00790423,1.45921e-06,-5.92451e-10,4.81983e-14,17936.6,-81.1625], Tmin=(1615.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]OC=C1CC([O])C1(22761)',
    structure = SMILES('[O]OC=C1CC([O])C1'),
    E0 = (241.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96412,0.0572769,-3.17668e-05,-1.42942e-10,4.13838e-12,29176.6,27.9907], Tmin=(100,'K'), Tmax=(1051.86,'K')), NASAPolynomial(coeffs=[14.4049,0.0239026,-9.46919e-06,1.75724e-09,-1.23745e-13,25367.7,-42.1975], Tmin=(1051.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([O])CC1=CO1(22841)',
    structure = SMILES('[CH2]C([O])CC1=CO1'),
    E0 = (336.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650842,0.0640753,-4.84619e-05,8.5732e-09,3.60294e-12,40658.2,26.9519], Tmin=(100,'K'), Tmax=(975.758,'K')), NASAPolynomial(coeffs=[17.7507,0.014908,-5.05628e-06,9.01901e-10,-6.42947e-14,36324.7,-60.2317], Tmin=(975.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C[C]=CO[O](22842)',
    structure = SMILES('C=C(O)C[C]=CO[O]'),
    E0 = (179.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0907145,0.0794693,-8.40538e-05,4.10118e-08,-6.66537e-12,21764.3,32.0462], Tmin=(100,'K'), Tmax=(959.79,'K')), NASAPolynomial(coeffs=[18.8239,0.0150679,-4.77003e-06,7.81932e-10,-5.20598e-14,17538.7,-60.8423], Tmin=(959.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(=O)C[C]=CO[O](22843)',
    structure = SMILES('CC(=O)C[C]=CO[O]'),
    E0 = (157.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3515,0.0623557,-5.23623e-05,2.2599e-08,-4.10263e-12,19039.2,27.8976], Tmin=(100,'K'), Tmax=(1255.38,'K')), NASAPolynomial(coeffs=[10.9864,0.031656,-1.56803e-05,3.11892e-09,-2.23287e-13,16620.2,-20.7764], Tmin=(1255.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])CC=CO[O](3038)',
    structure = SMILES('C=C([O])CC=CO[O]'),
    E0 = (79.7029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,236.266,236.266,236.266],'cm^-1')),
        HinderedRotor(inertia=(0.22083,'amu*angstrom^2'), symmetry=1, barrier=(8.74759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22083,'amu*angstrom^2'), symmetry=1, barrier=(8.74759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22083,'amu*angstrom^2'), symmetry=1, barrier=(8.74759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459127,0.0712632,-6.80877e-05,3.05547e-08,-4.4767e-12,9719.38,31.7509], Tmin=(100,'K'), Tmax=(992.009,'K')), NASAPolynomial(coeffs=[16.4305,0.0186644,-6.39843e-06,1.08943e-09,-7.32521e-14,5969.95,-48.1004], Tmin=(992.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.7029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(O)C=C=CO[O](19318)',
    structure = SMILES('[CH2]C(O)C=C=CO[O]'),
    E0 = (230.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230512,0.0815321,-9.54012e-05,5.68697e-08,-1.33029e-11,27865.3,33.0319], Tmin=(100,'K'), Tmax=(1048.36,'K')), NASAPolynomial(coeffs=[16.172,0.0207067,-8.37034e-06,1.5246e-09,-1.04691e-13,24522.9,-44.629], Tmin=(1048.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C=C=CO[O](19326)',
    structure = SMILES('CC([O])C=C=CO[O]'),
    E0 = (249.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.611814,'amu*angstrom^2'), symmetry=1, barrier=(14.0668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611895,'amu*angstrom^2'), symmetry=1, barrier=(14.0687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610694,'amu*angstrom^2'), symmetry=1, barrier=(14.0411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448872,0.0735047,-7.38026e-05,3.78601e-08,-7.68324e-12,30117.8,31.4207], Tmin=(100,'K'), Tmax=(1197.63,'K')), NASAPolynomial(coeffs=[16.0341,0.021452,-8.60931e-06,1.57063e-09,-1.08113e-13,26384.6,-46.5802], Tmin=(1197.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([O])C=C=COO(20483)',
    structure = SMILES('[CH2]C([O])C=C=COO'),
    E0 = (308.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0373882,0.0822197,-8.8565e-05,4.72656e-08,-9.85337e-12,37299.3,32.773], Tmin=(100,'K'), Tmax=(1173.99,'K')), NASAPolynomial(coeffs=[18.5583,0.0191159,-7.93779e-06,1.48047e-09,-1.03492e-13,32950.6,-59.5504], Tmin=(1173.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = '[CH2]C([O])C=C=CO[O](20484)',
    structure = SMILES('[CH2]C([O])C=C=CO[O]'),
    E0 = (460.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.598323,'amu*angstrom^2'), symmetry=1, barrier=(13.7566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.597998,'amu*angstrom^2'), symmetry=1, barrier=(13.7492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598081,'amu*angstrom^2'), symmetry=1, barrier=(13.751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381666,0.0776828,-9.00261e-05,5.25466e-08,-1.2006e-11,55566.3,32.2969], Tmin=(100,'K'), Tmax=(1073.44,'K')), NASAPolynomial(coeffs=[16.1915,0.0187694,-7.70088e-06,1.41731e-09,-9.80591e-14,52172.2,-45.0964], Tmin=(1073.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CC#CO[O](22844)',
    structure = SMILES('[CH2]C([O])CC#CO[O]'),
    E0 = (522.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2100,2250,500,550,248.768,248.875,2880.67],'cm^-1')),
        HinderedRotor(inertia=(0.339853,'amu*angstrom^2'), symmetry=1, barrier=(14.9231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338781,'amu*angstrom^2'), symmetry=1, barrier=(14.9185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67283,'amu*angstrom^2'), symmetry=1, barrier=(73.4239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67228,'amu*angstrom^2'), symmetry=1, barrier=(73.4513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461937,0.0805996,-0.000108138,7.82595e-08,-2.25143e-11,62954.9,31.9357], Tmin=(100,'K'), Tmax=(853.153,'K')), NASAPolynomial(coeffs=[12.3211,0.0249941,-1.03669e-05,1.85445e-09,-1.23801e-13,60931.5,-23.3938], Tmin=(853.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO)"""),
)

species(
    label = '[O]OC=[C]CC=O(21077)',
    structure = SMILES('[O]OC=[C]CC=O'),
    E0 = (212.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,1685,370,320.754],'cm^-1')),
        HinderedRotor(inertia=(0.163811,'amu*angstrom^2'), symmetry=1, barrier=(11.7042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163135,'amu*angstrom^2'), symmetry=1, barrier=(11.6947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(11.7126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93173,0.0479998,-3.92548e-05,1.57778e-08,-2.61749e-12,25619,24.2192], Tmin=(100,'K'), Tmax=(1373.67,'K')), NASAPolynomial(coeffs=[10.6644,0.0225711,-1.14876e-05,2.30183e-09,-1.64942e-13,23219.8,-20.6832], Tmin=(1373.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=CC[C]=CO[O](21013)',
    structure = SMILES('C=CC[C]=CO[O]'),
    E0 = (394.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.44098,'amu*angstrom^2'), symmetry=1, barrier=(10.139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440859,'amu*angstrom^2'), symmetry=1, barrier=(10.1362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440853,'amu*angstrom^2'), symmetry=1, barrier=(10.1361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23463,0.0590259,-5.532e-05,2.78936e-08,-5.71831e-12,47524.2,27.8847], Tmin=(100,'K'), Tmax=(1168.22,'K')), NASAPolynomial(coeffs=[11.6643,0.0233142,-9.46559e-06,1.72571e-09,-1.18294e-13,45087.4,-24.054], Tmin=(1168.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S)"""),
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
    label = '[CH]=[C]CC([CH2])[O](15961)',
    structure = SMILES('[CH]=[C]CC([CH2])[O]'),
    E0 = (714.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,358.675,359.013],'cm^-1')),
        HinderedRotor(inertia=(0.00130894,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134214,'amu*angstrom^2'), symmetry=1, barrier=(12.2609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134105,'amu*angstrom^2'), symmetry=1, barrier=(12.2629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.28011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34582,0.059919,-6.53842e-05,3.89315e-08,-9.43846e-12,85982.6,25.5326], Tmin=(100,'K'), Tmax=(994.722,'K')), NASAPolynomial(coeffs=[10.4569,0.0232814,-1.01361e-05,1.90406e-09,-1.325e-13,84170,-18.3748], Tmin=(994.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]([O])C[C]=CO[O](22845)',
    structure = SMILES('[CH2][C]([O])C[C]=CO[O]'),
    E0 = (711.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,228.718,229.395,230.752],'cm^-1')),
        HinderedRotor(inertia=(0.170273,'amu*angstrom^2'), symmetry=1, barrier=(6.39325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174618,'amu*angstrom^2'), symmetry=1, barrier=(6.38517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171358,'amu*angstrom^2'), symmetry=1, barrier=(6.37779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170339,'amu*angstrom^2'), symmetry=1, barrier=(6.38711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00638986,0.0951765,-0.000153647,1.27629e-07,-4.09497e-11,85693.7,34.0617], Tmin=(100,'K'), Tmax=(873.315,'K')), NASAPolynomial(coeffs=[12.3861,0.0254242,-1.15262e-05,2.10322e-09,-1.39743e-13,84024.6,-21.2096], Tmin=(873.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(711.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C=[C][CH]O[O](20487)',
    structure = SMILES('[CH2]C([O])C=[C][CH]O[O]'),
    E0 = (642.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,180,896.23,896.847],'cm^-1')),
        HinderedRotor(inertia=(0.131879,'amu*angstrom^2'), symmetry=1, barrier=(3.03217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708631,'amu*angstrom^2'), symmetry=1, barrier=(16.2928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0429653,'amu*angstrom^2'), symmetry=1, barrier=(24.6619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709308,'amu*angstrom^2'), symmetry=1, barrier=(16.3084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580765,0.0724438,-7.55185e-05,3.97289e-08,-8.26929e-12,77353.1,34.214], Tmin=(100,'K'), Tmax=(1166.65,'K')), NASAPolynomial(coeffs=[15.6795,0.0206759,-8.95899e-06,1.69436e-09,-1.18913e-13,73830.1,-40.9559], Tmin=(1166.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(642.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C[C]=[C]O[O](22846)',
    structure = SMILES('[CH2]C([O])C[C]=[C]O[O]'),
    E0 = (774.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,180.033,2443.71,2444.49],'cm^-1')),
        HinderedRotor(inertia=(0.505921,'amu*angstrom^2'), symmetry=1, barrier=(12.1194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502606,'amu*angstrom^2'), symmetry=1, barrier=(12.1123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.505149,'amu*angstrom^2'), symmetry=1, barrier=(12.1092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504612,'amu*angstrom^2'), symmetry=1, barrier=(12.1155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367106,0.088094,-0.000143149,1.23859e-07,-4.14908e-11,93270.2,36.425], Tmin=(100,'K'), Tmax=(863.86,'K')), NASAPolynomial(coeffs=[9.53273,0.0296274,-1.38002e-05,2.5608e-09,-1.72249e-13,92284.7,-2.99132], Tmin=(863.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(774.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](O)C[C]=CO[O](22847)',
    structure = SMILES('[CH2][C](O)C[C]=CO[O]'),
    E0 = (480.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.170056,0.0991787,-0.000159569,1.32653e-07,-4.25304e-11,57993.2,34.841], Tmin=(100,'K'), Tmax=(880.227,'K')), NASAPolynomial(coeffs=[12.4742,0.027177,-1.20885e-05,2.18507e-09,-1.44256e-13,56330.6,-21.3471], Tmin=(880.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C[C]([O])C[C]=CO[O](22848)',
    structure = SMILES('C[C]([O])C[C]=CO[O]'),
    E0 = (499.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.183517,0.0895297,-0.000132249,1.0631e-07,-3.38874e-11,60239.9,32.7474], Tmin=(100,'K'), Tmax=(838.515,'K')), NASAPolynomial(coeffs=[11.4071,0.029493,-1.32299e-05,2.44374e-09,-1.65284e-13,58586.1,-18.0614], Tmin=(838.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C=C[CH]O[O](20490)',
    structure = SMILES('[CH2]C([O])C=C[CH]O[O]'),
    E0 = (404.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.191959,0.0726724,-6.6926e-05,3.0473e-08,-5.42037e-12,48768.2,35.1447], Tmin=(100,'K'), Tmax=(1370.84,'K')), NASAPolynomial(coeffs=[19.0478,0.0176528,-6.72282e-06,1.19513e-09,-8.10067e-14,43598.5,-61.7716], Tmin=(1370.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CC=[C]O[O](22849)',
    structure = SMILES('[CH2]C([O])CC=[C]O[O]'),
    E0 = (536.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461501,0.0827192,-0.000115414,9.0442e-08,-2.8623e-11,64664.3,35.616], Tmin=(100,'K'), Tmax=(800.669,'K')), NASAPolynomial(coeffs=[10.492,0.0305349,-1.37645e-05,2.57015e-09,-1.75838e-13,63124.6,-10.1302], Tmin=(800.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([O])C[C]=[C]OO(22850)',
    structure = SMILES('[CH2]C([O])C[C]=[C]OO'),
    E0 = (622.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146104,0.0911572,-0.000136527,1.12068e-07,-3.6732e-11,74998,36.4612], Tmin=(100,'K'), Tmax=(809.634,'K')), NASAPolynomial(coeffs=[11.0231,0.0314626,-1.48955e-05,2.82675e-09,-1.94507e-13,73432,-12.5115], Tmin=(809.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)C=[C][CH]O[O](20496)',
    structure = SMILES('[CH2]C(O)C=[C][CH]O[O]'),
    E0 = (411.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443871,0.0761276,-8.03334e-05,4.33624e-08,-9.29103e-12,49651.5,34.8978], Tmin=(100,'K'), Tmax=(1133.25,'K')), NASAPolynomial(coeffs=[15.5573,0.0227818,-9.72311e-06,1.82358e-09,-1.27337e-13,46226,-39.9063], Tmin=(1133.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(CJCO) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])C=[C][CH]O[O](20497)',
    structure = SMILES('CC([O])C=[C][CH]O[O]'),
    E0 = (430.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,332.938,332.942,332.945],'cm^-1')),
        HinderedRotor(inertia=(0.00152085,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00152087,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.277104,'amu*angstrom^2'), symmetry=1, barrier=(21.7959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.277103,'amu*angstrom^2'), symmetry=1, barrier=(21.7958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574739,0.0690743,-6.18582e-05,2.79868e-08,-5.04266e-12,51907.8,33.604], Tmin=(100,'K'), Tmax=(1333.4,'K')), NASAPolynomial(coeffs=[16.1386,0.0223854,-9.33623e-06,1.72729e-09,-1.19293e-13,47757.2,-45.9609], Tmin=(1333.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C[CH][CH]O[O](3040)',
    structure = SMILES('C=C([O])C[CH][CH]O[O]'),
    E0 = (335.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,342.149,343.271,344.524,1332.43],'cm^-1')),
        HinderedRotor(inertia=(0.0939153,'amu*angstrom^2'), symmetry=1, barrier=(7.85382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0941713,'amu*angstrom^2'), symmetry=1, barrier=(7.87926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0938533,'amu*angstrom^2'), symmetry=1, barrier=(7.86464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00142827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263024,0.0876436,-0.000127135,1.01888e-07,-3.2548e-11,40495.6,33.9741], Tmin=(100,'K'), Tmax=(835.222,'K')), NASAPolynomial(coeffs=[10.8878,0.0306165,-1.36849e-05,2.52699e-09,-1.71067e-13,38935.1,-14.0886], Tmin=(835.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C(O)C[C]=[C]O[O](22851)',
    structure = SMILES('[CH2]C(O)C[C]=[C]O[O]'),
    E0 = (544.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.204892,0.0920768,-0.000148991,1.28758e-07,-4.30059e-11,65569.7,37.1993], Tmin=(100,'K'), Tmax=(870.405,'K')), NASAPolynomial(coeffs=[9.62099,0.0313802,-1.43626e-05,2.64269e-09,-1.76766e-13,64590.6,-3.12938], Tmin=(870.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CJCO) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'CC([O])C[C]=[C]O[O](22852)',
    structure = SMILES('CC([O])C[C]=[C]O[O]'),
    E0 = (562.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,287.834,287.836,3568.96],'cm^-1')),
        HinderedRotor(inertia=(0.157871,'amu*angstrom^2'), symmetry=1, barrier=(9.28171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157866,'amu*angstrom^2'), symmetry=1, barrier=(9.28168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157866,'amu*angstrom^2'), symmetry=1, barrier=(9.28178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157875,'amu*angstrom^2'), symmetry=1, barrier=(9.28174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.551363,0.0825224,-0.000122062,1.03026e-07,-3.4682e-11,67816.7,35.1305], Tmin=(100,'K'), Tmax=(834.859,'K')), NASAPolynomial(coeffs=[8.55636,0.0336911,-1.55006e-05,2.9005e-09,-1.97719e-13,66845.3,0.142719], Tmin=(834.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([O])C=[C][CH]OO(20492)',
    structure = SMILES('[CH2]C([O])C=[C][CH]OO'),
    E0 = (490.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.180552,0.0776012,-7.60381e-05,3.67434e-08,-6.97988e-12,59088.5,34.8933], Tmin=(100,'K'), Tmax=(1277.6,'K')), NASAPolynomial(coeffs=[18.4631,0.0203609,-8.83339e-06,1.67511e-09,-1.1771e-13,54417,-57.788], Tmin=(1277.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])C[C]=COO(22853)',
    structure = SMILES('[CH2][C]([O])C[C]=COO'),
    E0 = (559.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0551857,0.0958601,-0.000136731,9.88125e-08,-2.67606e-11,67414.1,33.5004], Tmin=(100,'K'), Tmax=(694.865,'K')), NASAPolynomial(coeffs=[13.684,0.0276202,-1.28441e-05,2.42437e-09,-1.66746e-13,65242.8,-29.6666], Tmin=(694.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S)"""),
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
    E0 = (534.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (692.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (704.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (534.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (642.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1184.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (973.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1074.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (983.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (540.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (543.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (543.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (543.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (666.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (557.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (557.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (623.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (612.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (612.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (559.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (550.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (681.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (749.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (650.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (637.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (534.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (635.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (534.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (705.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (805.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (923.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (853.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (986.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (648.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (676.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (692.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (768.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (759.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (610.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (652.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (676.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (659.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (595.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (585.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (633.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['vinoxy(1351)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[O][CH]CC[C]=CO[O](21034)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['O2(S)(666)', 'C#CCC([CH2])[O](15960)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2]C([O])C[C]=C[O](22838)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2][CH]C[C]=CO[O](21352)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[O][CH]C[C]=CO[O](22760)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]=CO[O](21346)', '[CH2]C([CH2])[O](1366)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([O])C[C]=CO[O](21141)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[O]OC=[C]CC1CO1(22839)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])CC1=COO1(22840)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C1CC(=CO[O])O1(22821)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[O]OC=C1CC([O])C1(22761)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['O(4)', '[CH2]C([O])CC1=CO1(22841)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO;Y_rad_intra;OOJ] for rate rule [R2OO_D;Cd_rad_in/NonDeC;OOJ]
Euclidian distance = 3.1622776601683795
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['C=C(O)C[C]=CO[O](22842)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['CC(=O)C[C]=CO[O](22843)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['C=C([O])CC=CO[O](3038)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C(O)C=C=CO[O](19318)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['CC([O])C=C=CO[O](19326)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])C=C=COO(20483)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', 'C=C([O])C[C]=CO[O](21045)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C([O])C=C=CO[O](20484)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C([O])CC#CO[O](22844)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(20)', '[O]OC=[C]CC=O(21077)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)', 'C=CC[C]=CO[O](21013)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['vinoxy(1351)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(90.8835,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 87.2 to 90.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH][O](1556)', 'C=C=CO[O](16806)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O2(2)', 'C#CCC([CH2])[O](15960)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.3016e+08,'cm^3/(mol*s)'), n=1.628, Ea=(148.067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-Cs;OJ] for rate rule [Ct-H_Ct-Cs;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.9 to 148.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['O2(2)', '[CH]=[C]CC([CH2])[O](15961)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
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
    reactants = ['[CH2][CH][O](1556)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][C]([O])C[C]=CO[O](22845)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]C([O])C=[C][CH]O[O](20487)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2]C([O])C[C]=[C]O[O](22846)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2][C](O)C[C]=CO[O](22847)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['C[C]([O])C[C]=CO[O](22848)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])C=C[CH]O[O](20490)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])CC=[C]O[O](22849)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([O])C[C]=[C]OO(22850)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.58963e+06,'s^-1'), n=1.78157, Ea=(137.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C(O)C=[C][CH]O[O](20496)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['CC([O])C=[C][CH]O[O](20497)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['C=C([O])C[CH][CH]O[O](3040)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C(O)C[C]=[C]O[O](22851)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(49.9162,'s^-1'), n=2.87175, Ea=(124.402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cd_H_out_single] for rate rule [R5HJ_3;O_rad_out;Cd_H_out_singleNd]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CC([O])C[C]=[C]O[O](22852)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])C=[C][CH]OO(20492)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C]([O])C[C]=COO(22853)'],
    products = ['[CH2]C([O])C[C]=CO[O](21021)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4696',
    isomers = [
        '[CH2]C([O])C[C]=CO[O](21021)',
    ],
    reactants = [
        ('vinoxy(1351)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4696',
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


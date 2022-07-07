species(
    label = 'C=C([CH]O[O])C[CH][O](20382)',
    structure = SMILES('C=C([CH]O[O])C[CH][O]'),
    E0 = (378.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,263.452,263.536,264.358,1655.34],'cm^-1')),
        HinderedRotor(inertia=(0.00688885,'amu*angstrom^2'), symmetry=1, barrier=(45.1789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917495,'amu*angstrom^2'), symmetry=1, barrier=(45.1777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275991,'amu*angstrom^2'), symmetry=1, barrier=(13.6879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799187,0.0776825,-9.30708e-05,5.78414e-08,-1.12702e-11,45603.6,31.3517], Tmin=(100,'K'), Tmax=(620.58,'K')), NASAPolynomial(coeffs=[8.92717,0.0357517,-1.69998e-05,3.27829e-09,-2.29702e-13,44393.4,-5.60588], Tmin=(620.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJO) + radical(CCsJOH)"""),
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
    label = 'C=C([CH]O[O])CC=O(21027)',
    structure = SMILES('C=C([CH]O[O])CC=O'),
    E0 = (55.4753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812471,0.0596919,-3.59132e-05,5.36606e-09,1.2423e-12,6795.45,29.7553], Tmin=(100,'K'), Tmax=(1237.44,'K')), NASAPolynomial(coeffs=[16.5438,0.0242166,-1.15491e-05,2.28131e-09,-1.63142e-13,1724.92,-54.2473], Tmin=(1237.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.4753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO)"""),
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
    label = 'C=[C]C(C[CH][O])O[O](21033)',
    structure = SMILES('C=[C]C(C[CH][O])O[O]'),
    E0 = (494.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,298.46,298.46,298.46,1468.41],'cm^-1')),
        HinderedRotor(inertia=(0.115145,'amu*angstrom^2'), symmetry=1, barrier=(7.2786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115146,'amu*angstrom^2'), symmetry=1, barrier=(7.2786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115146,'amu*angstrom^2'), symmetry=1, barrier=(7.27859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115144,'amu*angstrom^2'), symmetry=1, barrier=(7.27858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4357.22,'J/mol'), sigma=(7.17033,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.59 K, Pc=26.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272119,0.0928565,-0.00015423,1.41193e-07,-5.01045e-11,59588.6,34.0765], Tmin=(100,'K'), Tmax=(844.441,'K')), NASAPolynomial(coeffs=[7.04941,0.0387588,-1.90652e-05,3.63926e-09,-2.49475e-13,59228.2,7.16919], Tmin=(844.441,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH) + radical(Cds_S)"""),
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
    label = 'C=C([CH][O])C[CH][O](22742)',
    structure = SMILES('[CH2]C(=C[O])C[CH][O]'),
    E0 = (274.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,412.819,412.824,412.827,412.839],'cm^-1')),
        HinderedRotor(inertia=(0.000989139,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150504,'amu*angstrom^2'), symmetry=1, barrier=(18.2017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353698,'amu*angstrom^2'), symmetry=1, barrier=(42.7746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736661,0.070042,-7.29009e-05,3.95087e-08,-8.54868e-12,33166,27.2138], Tmin=(100,'K'), Tmax=(1120.14,'K')), NASAPolynomial(coeffs=[14.0047,0.0226622,-9.45375e-06,1.7473e-09,-1.2087e-13,30193.5,-38.3022], Tmin=(1120.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CCsJOH) + radical(Allyl_P)"""),
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
    label = '[CH]C(=C)C[CH][O](15924)',
    structure = SMILES('[CH]C(=C)C[CH][O]'),
    E0 = (561.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,538.002,538.002,538.003,538.003,538.003,538.003],'cm^-1')),
        HinderedRotor(inertia=(0.261593,'amu*angstrom^2'), symmetry=1, barrier=(53.7307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261593,'amu*angstrom^2'), symmetry=1, barrier=(53.7307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261593,'amu*angstrom^2'), symmetry=1, barrier=(53.7307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.1,'J/mol'), sigma=(6.27192,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.76 K, Pc=33.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35295,0.0638251,-8.08714e-05,6.93203e-08,-2.50285e-11,67596.7,24.5591], Tmin=(100,'K'), Tmax=(805.75,'K')), NASAPolynomial(coeffs=[4.47181,0.0392434,-1.81715e-05,3.42873e-09,-2.36098e-13,67389.4,12.0189], Tmin=(805.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C(=C)[CH]O[O](21386)',
    structure = SMILES('[CH2]C(=C)[CH]O[O]'),
    E0 = (306.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,983.224],'cm^-1')),
        HinderedRotor(inertia=(0.266225,'amu*angstrom^2'), symmetry=1, barrier=(21.8329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22972,'amu*angstrom^2'), symmetry=1, barrier=(28.2736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118064,'amu*angstrom^2'), symmetry=1, barrier=(80.835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6098,0.0438254,-1.89927e-05,-9.76165e-09,7.61582e-12,36925.6,22.6552], Tmin=(100,'K'), Tmax=(991.431,'K')), NASAPolynomial(coeffs=[13.3403,0.015827,-5.87653e-06,1.0834e-09,-7.75518e-14,33649.6,-38.6276], Tmin=(991.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
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
    label = 'C=[C]C[CH][O](2773)',
    structure = SMILES('C=[C]C[CH][O]'),
    E0 = (467.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,180,180,1699.6],'cm^-1')),
        HinderedRotor(inertia=(0.204683,'amu*angstrom^2'), symmetry=1, barrier=(4.70606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205128,'amu*angstrom^2'), symmetry=1, barrier=(4.71629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07829,0.0495382,-7.92644e-05,7.71072e-08,-2.8994e-11,56288.9,20.5615], Tmin=(100,'K'), Tmax=(845.873,'K')), NASAPolynomial(coeffs=[2.88158,0.0292521,-1.4053e-05,2.66818e-09,-1.82773e-13,56742.8,20.3072], Tmin=(845.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S)"""),
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
    label = 'C=C([C]O[O])C[CH][O](22743)',
    structure = SMILES('[CH2]C(=[C]O[O])C[CH][O]'),
    E0 = (649.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,180,1430.38,1433.24],'cm^-1')),
        HinderedRotor(inertia=(0.288212,'amu*angstrom^2'), symmetry=1, barrier=(6.62655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286748,'amu*angstrom^2'), symmetry=1, barrier=(6.59291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287303,'amu*angstrom^2'), symmetry=1, barrier=(6.60565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287112,'amu*angstrom^2'), symmetry=1, barrier=(6.60128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39532,0.0918963,-0.00016239,1.51701e-07,-5.36078e-11,78240.4,34.9559], Tmin=(100,'K'), Tmax=(871.075,'K')), NASAPolynomial(coeffs=[6.26824,0.0363851,-1.7649e-05,3.3083e-09,-2.22785e-13,78300.2,13.6487], Tmin=(871.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'C=C([CH]O[O])C[C][O](22744)',
    structure = SMILES('C=C([CH]O[O])C[C][O]'),
    E0 = (658.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769401,0.0763574,-9.49771e-05,6.51254e-08,-1.83263e-11,79366.1,30.6171], Tmin=(100,'K'), Tmax=(858.238,'K')), NASAPolynomial(coeffs=[10.5998,0.0305413,-1.49019e-05,2.92478e-09,-2.07706e-13,77678.7,-15.306], Tmin=(858.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(658.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJO) + radical(CH2_triplet)"""),
)

species(
    label = 'C=C(C[CH][O])C1OO1(22745)',
    structure = SMILES('C=C(C[CH][O])C1OO1'),
    E0 = (224.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.653968,0.0727478,-7.56798e-05,4.27255e-08,-9.74173e-12,27152.6,29.3726], Tmin=(100,'K'), Tmax=(1061.1,'K')), NASAPolynomial(coeffs=[12.836,0.0268247,-1.07605e-05,1.93726e-09,-1.31654e-13,24567.4,-30.1208], Tmin=(1061.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(dioxirane) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=C1CC([O])C1O[O](22746)',
    structure = SMILES('C=C1CC([O])C1O[O]'),
    E0 = (224.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930326,0.0575014,-3.11058e-05,-1.39047e-10,3.70611e-12,27133,26.8562], Tmin=(100,'K'), Tmax=(1094.92,'K')), NASAPolynomial(coeffs=[14.7451,0.0246985,-1.03681e-05,1.96957e-09,-1.3979e-13,23048.9,-45.8803], Tmin=(1094.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C(C=C[O])CO[O](10577)',
    structure = SMILES('C=C(C=C[O])CO[O]'),
    E0 = (59.3745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.24509,0.0793835,-8.06162e-05,4.02687e-08,-7.69029e-12,7305.84,30.383], Tmin=(100,'K'), Tmax=(1385.45,'K')), NASAPolynomial(coeffs=[20.8376,0.0133608,-3.55481e-06,5.024e-10,-3.00849e-14,1958.68,-76.417], Tmin=(1385.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.3745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([CH][CH][O])[CH]O[O](22747)',
    structure = SMILES('[CH2]C([CH][CH][O])[CH]O[O]'),
    E0 = (729.646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,201.802,801.2,1205.4,1609.6],'cm^-1')),
        HinderedRotor(inertia=(0.154002,'amu*angstrom^2'), symmetry=1, barrier=(3.58341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154002,'amu*angstrom^2'), symmetry=1, barrier=(3.58341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154002,'amu*angstrom^2'), symmetry=1, barrier=(3.58341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154002,'amu*angstrom^2'), symmetry=1, barrier=(3.58341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154002,'amu*angstrom^2'), symmetry=1, barrier=(3.58341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.161007,0.0982959,-0.000173223,1.61358e-07,-5.64895e-11,87881.2,37.1128], Tmin=(100,'K'), Tmax=(886.636,'K')), NASAPolynomial(coeffs=[5.96533,0.0394384,-1.83753e-05,3.36643e-09,-2.22969e-13,88136.1,17.0506], Tmin=(886.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(ROOJ) + radical(CCJCO) + radical(CCsJOOH) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C[C]1CC1O[O](22748)',
    structure = SMILES('[O][CH]C[C]1CC1O[O]'),
    E0 = (473.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5859,0.0572249,-4.57077e-05,2.08523e-08,-4.22784e-12,57043.9,31.632], Tmin=(100,'K'), Tmax=(1104.48,'K')), NASAPolynomial(coeffs=[7.54667,0.0356375,-1.63902e-05,3.15646e-09,-2.2242e-13,55727.2,2.28232], Tmin=(1104.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopropane) + radical(CCOJ) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOH)"""),
)

species(
    label = '[O]O[CH][C]1CC([O])C1(22749)',
    structure = SMILES('[O]O[CH][C]1CC([O])C1'),
    E0 = (482.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37879,0.0477046,-1.06669e-05,-1.61602e-08,8.26863e-12,58088.1,31.9331], Tmin=(100,'K'), Tmax=(1071.67,'K')), NASAPolynomial(coeffs=[12.1206,0.0280023,-1.16314e-05,2.195e-09,-1.55263e-13,54614.9,-26.0965], Tmin=(1071.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOOH)"""),
)

species(
    label = '[O]O[CH][C]1C[CH]OC1(22750)',
    structure = SMILES('[O]O[CH][C]1C[CH]OC1'),
    E0 = (392.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43593,0.0490215,-9.95643e-06,-2.6982e-08,1.68975e-11,47275.9,28.0468], Tmin=(100,'K'), Tmax=(863.917,'K')), NASAPolynomial(coeffs=[11.815,0.0239068,-6.18301e-06,8.44156e-10,-4.97824e-14,44626.5,-25.463], Tmin=(863.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOCs) + radical(CCsJOOH)"""),
)

species(
    label = '[O][CH]C[C]1[CH]OOC1(22751)',
    structure = SMILES('[O][CH]C[C]1[CH]OOC1'),
    E0 = (442.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37229,0.0546049,-3.374e-05,5.94629e-09,2.47018e-12,53318.6,30.8484], Tmin=(100,'K'), Tmax=(883.848,'K')), NASAPolynomial(coeffs=[9.02122,0.0304244,-1.04137e-05,1.71085e-09,-1.10479e-13,51558.9,-7.41476], Tmin=(883.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(12dioxolane) + radical(CCOJ) + radical(C2CJCOOH) + radical(CCsJOOC) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([CH]O[O])CC1[O](22752)',
    structure = SMILES('[CH2]C1([CH]O[O])CC1[O]'),
    E0 = (488.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.417308,0.0717682,-6.588e-05,3.04235e-08,-5.56516e-12,58895.5,29.5457], Tmin=(100,'K'), Tmax=(1320.05,'K')), NASAPolynomial(coeffs=[16.9733,0.0216002,-8.87304e-06,1.63317e-09,-1.12641e-13,54524.5,-54.9247], Tmin=(1320.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(ROOJ) + radical(CCsJOOH) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C1([CH]O[O])C[CH]O1(22753)',
    structure = SMILES('[CH2]C1([CH]O[O])C[CH]O1'),
    E0 = (468.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.40306,0.0955639,-0.000131276,9.13342e-08,-2.38807e-11,56534.2,28.7484], Tmin=(100,'K'), Tmax=(1089.38,'K')), NASAPolynomial(coeffs=[16.4583,0.0182551,-3.62692e-06,2.42528e-10,1.18671e-15,53774.2,-49.848], Tmin=(1089.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(CCsJOOH) + radical(CCsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1(C[CH][O])[CH]OO1(22754)',
    structure = SMILES('[CH2]C1(C[CH][O])[CH]OO1'),
    E0 = (529.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.182416,0.0948972,-0.000156056,1.41966e-07,-5.00565e-11,63860,28.2687], Tmin=(100,'K'), Tmax=(850.163,'K')), NASAPolynomial(coeffs=[7.11434,0.0398306,-1.92852e-05,3.65178e-09,-2.49076e-13,63492.8,0.723297], Tmin=(850.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(CCOJ) + radical(CCsJOO) + radical(CJCOOH) + radical(CCsJOH)"""),
)

species(
    label = 'C=C([CH]O[O])C=C[O](21133)',
    structure = SMILES('C=C([CH]O[O])C=C[O]'),
    E0 = (176.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24169,'amu*angstrom^2'), symmetry=1, barrier=(28.5489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24153,'amu*angstrom^2'), symmetry=1, barrier=(28.5453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24162,'amu*angstrom^2'), symmetry=1, barrier=(28.5472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.365176,0.063671,-2.47548e-05,-3.03472e-08,2.09388e-11,21394.3,28.7847], Tmin=(100,'K'), Tmax=(932.088,'K')), NASAPolynomial(coeffs=[23.1982,0.00714195,-5.00442e-07,2.40051e-11,-6.15085e-15,15336.9,-89.4263], Tmin=(932.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C=C(C=O)C[CH][O](22755)',
    structure = SMILES('C=C(C=O)C[CH][O]'),
    E0 = (98.8479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,180,180,710.594],'cm^-1')),
        HinderedRotor(inertia=(0.171464,'amu*angstrom^2'), symmetry=1, barrier=(3.9423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165921,'amu*angstrom^2'), symmetry=1, barrier=(3.81485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16327,'amu*angstrom^2'), symmetry=1, barrier=(3.75389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11087,0.0725699,-0.000113227,1.09218e-07,-4.22442e-11,11983.9,26.1813], Tmin=(100,'K'), Tmax=(789.364,'K')), NASAPolynomial(coeffs=[3.88533,0.0413952,-2.1463e-05,4.24949e-09,-2.99842e-13,12079.1,16.83], Tmin=(789.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.8479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = 'C=C([CH][CH][O])[CH]O[O](22756)',
    structure = SMILES('[CH2]C([CH]O[O])=C[CH][O]'),
    E0 = (460.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,429.522,431.684,432.491],'cm^-1')),
        HinderedRotor(inertia=(0.00091417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000921898,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000915935,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.371673,'amu*angstrom^2'), symmetry=1, barrier=(49.1652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71968,0.063026,-5.14009e-05,2.07409e-08,-3.32168e-12,55460,32.8839], Tmin=(100,'K'), Tmax=(1492.51,'K')), NASAPolynomial(coeffs=[16.684,0.0202408,-8.40085e-06,1.53377e-09,-1.04417e-13,50694.6,-50.5276], Tmin=(1492.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C([CH]O[O])C[CH][O](20379)',
    structure = SMILES('[CH]=C([CH]O[O])C[CH][O]'),
    E0 = (625.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3120,650,792.5,1650,322.897,323.665,2217.39],'cm^-1')),
        HinderedRotor(inertia=(0.00315042,'amu*angstrom^2'), symmetry=1, barrier=(10.9804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147662,'amu*angstrom^2'), symmetry=1, barrier=(10.987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56577,'amu*angstrom^2'), symmetry=1, barrier=(41.8008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564834,'amu*angstrom^2'), symmetry=1, barrier=(41.806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410394,0.0856238,-0.000129291,1.08654e-07,-3.65243e-11,75336.8,33.1095], Tmin=(100,'K'), Tmax=(810.2,'K')), NASAPolynomial(coeffs=[9.68124,0.0319788,-1.53948e-05,2.93943e-09,-2.02881e-13,74093,-8.0708], Tmin=(810.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJO) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = 'C=C([CH]C[O])[CH]O[O](22757)',
    structure = SMILES('C=C([CH]C[O])[CH]O[O]'),
    E0 = (314.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858031,0.0733087,-7.41073e-05,4.14787e-08,-9.75043e-12,37981.8,28.7628], Tmin=(100,'K'), Tmax=(1004.98,'K')), NASAPolynomial(coeffs=[10.5721,0.0346441,-1.63968e-05,3.19499e-09,-2.26766e-13,36029.3,-18.1503], Tmin=(1004.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJCO) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C]([CH]C=O)CO[O](10571)',
    structure = SMILES('[CH2][C](C=C[O])CO[O]'),
    E0 = (305.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105453,0.074288,-7.26244e-05,3.62388e-08,-7.02157e-12,36921.2,34.2221], Tmin=(100,'K'), Tmax=(1334.67,'K')), NASAPolynomial(coeffs=[18.0329,0.0174099,-5.16067e-06,7.72488e-10,-4.71082e-14,32416.3,-56.3917], Tmin=(1334.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(C[CH][O])CO[O](19320)',
    structure = SMILES('[CH]=C(C[CH][O])CO[O]'),
    E0 = (508.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,1734.9],'cm^-1')),
        HinderedRotor(inertia=(0.254577,'amu*angstrom^2'), symmetry=1, barrier=(5.85322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254617,'amu*angstrom^2'), symmetry=1, barrier=(5.85414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254514,'amu*angstrom^2'), symmetry=1, barrier=(5.85177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254578,'amu*angstrom^2'), symmetry=1, barrier=(5.85325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4356.36,'J/mol'), sigma=(7.16244,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.45 K, Pc=26.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.222468,0.0964823,-0.000168758,1.58797e-07,-5.67301e-11,61229.7,33.1842], Tmin=(100,'K'), Tmax=(864.451,'K')), NASAPolynomial(coeffs=[5.76654,0.0407596,-1.98922e-05,3.75263e-09,-2.542e-13,61394.7,13.7432], Tmin=(864.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = 'C=C([CH][CH]O)[CH]O[O](22758)',
    structure = SMILES('[CH2]C([CH]O[O])=C[CH]O'),
    E0 = (234.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521254,0.064091,-3.44886e-05,-8.35129e-09,9.56495e-12,28324.6,32.7202], Tmin=(100,'K'), Tmax=(984.776,'K')), NASAPolynomial(coeffs=[18.6121,0.0181992,-6.61225e-06,1.22724e-09,-8.91833e-14,23423.7,-61.0728], Tmin=(984.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])CC[O](20384)',
    structure = SMILES('[CH]=C([CH]O[O])CC[O]'),
    E0 = (445.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,350.032,350.293,2451.96],'cm^-1')),
        HinderedRotor(inertia=(0.128185,'amu*angstrom^2'), symmetry=1, barrier=(11.1585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00261544,'amu*angstrom^2'), symmetry=1, barrier=(11.1581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486997,'amu*angstrom^2'), symmetry=1, barrier=(42.4181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.487412,'amu*angstrom^2'), symmetry=1, barrier=(42.418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85327,0.0746099,-8.8995e-05,6.28002e-08,-1.86809e-11,53637.6,31.9856], Tmin=(100,'K'), Tmax=(806.348,'K')), NASAPolynomial(coeffs=[8.65266,0.03592,-1.70225e-05,3.29539e-09,-2.32082e-13,52379.7,-3.9632], Tmin=(806.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C[CH]O(20385)',
    structure = SMILES('[CH]=C([CH]O[O])C[CH]O'),
    E0 = (399.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3120,650,792.5,1650,180,680.455],'cm^-1')),
        HinderedRotor(inertia=(0.52302,'amu*angstrom^2'), symmetry=1, barrier=(12.0253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0364886,'amu*angstrom^2'), symmetry=1, barrier=(12.0267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928039,'amu*angstrom^2'), symmetry=1, barrier=(2.13374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523085,'amu*angstrom^2'), symmetry=1, barrier=(12.0268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237105,'amu*angstrom^2'), symmetry=1, barrier=(77.3628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.188037,0.0865927,-0.000110098,7.34456e-08,-1.94689e-11,48202.8,33.0553], Tmin=(100,'K'), Tmax=(922.73,'K')), NASAPolynomial(coeffs=[14.194,0.0258766,-1.13959e-05,2.13273e-09,-1.47496e-13,45618.1,-33.3887], Tmin=(922.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCsJOH) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = 'C=C([CH][CH][O])[CH]OO(22759)',
    structure = SMILES('[CH2]C([CH]OO)=C[CH][O]'),
    E0 = (308.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50659,0.0661085,-4.53191e-05,1.01794e-08,7.70764e-13,37187.2,32.8836], Tmin=(100,'K'), Tmax=(1135.4,'K')), NASAPolynomial(coeffs=[17.0425,0.0237517,-1.03652e-05,1.98884e-09,-1.41381e-13,32407.4,-53.5054], Tmin=(1135.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C([CH]OO)C[CH][O](20381)',
    structure = SMILES('[CH]=C([CH]OO)C[CH][O]'),
    E0 = (473.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3120,650,792.5,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.348185,0.086599,-0.000114251,8.41136e-08,-2.5384e-11,57057.8,32.5893], Tmin=(100,'K'), Tmax=(803.286,'K')), NASAPolynomial(coeffs=[10.8208,0.0344506,-1.68744e-05,3.29923e-09,-2.33117e-13,55375.3,-15.6414], Tmin=(803.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CCJO) + radical(CCsJOH) + radical(Cds_P)"""),
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
    label = '[O][CH]CC1=COOC1(19175)',
    structure = SMILES('[O][CH]CC1=COOC1'),
    E0 = (190.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21337,0.0618245,-5.0009e-05,2.1383e-08,-3.82648e-12,22966.5,28.2875], Tmin=(100,'K'), Tmax=(1287.43,'K')), NASAPolynomial(coeffs=[11.2223,0.0307273,-1.37773e-05,2.62124e-09,-1.83226e-13,20389.4,-22.5283], Tmin=(1287.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1=COOC([O])C1(22762)',
    structure = SMILES('C=C1[CH]OOC([O])C1'),
    E0 = (70.4789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50316,0.0337268,5.28795e-05,-8.77485e-08,3.29809e-11,8584.84,23.885], Tmin=(100,'K'), Tmax=(1045.14,'K')), NASAPolynomial(coeffs=[16.804,0.0277176,-1.39179e-05,2.96979e-09,-2.27218e-13,2516.44,-64.3391], Tmin=(1045.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.4789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(CCOJ) + radical(C=CCJO)"""),
)

species(
    label = '[O][CH]CC1=COC1(20375)',
    structure = SMILES('[O][CH]CC1=COC1'),
    E0 = (189.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982732,0.0546101,-2.22717e-05,-1.67379e-08,1.20656e-11,22886.1,23.6575], Tmin=(100,'K'), Tmax=(968.351,'K')), NASAPolynomial(coeffs=[16.8969,0.016329,-5.50335e-06,9.97984e-10,-7.25723e-14,18516.7,-59.2541], Tmin=(968.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1=COC([O])C1(22763)',
    structure = SMILES('[CH2]C1=COC([O])C1'),
    E0 = (9.39456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71081,0.0314799,4.70084e-05,-9.29195e-08,4.10593e-11,1229.85,18.9208], Tmin=(100,'K'), Tmax=(918.895,'K')), NASAPolynomial(coeffs=[17.3767,0.0121996,-1.36549e-06,1.0594e-10,-1.03838e-14,-3714.28,-66.5696], Tmin=(918.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.39456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = 'CC(C=C[O])=CO[O](22764)',
    structure = SMILES('CC([CH]O[O])=CC=O'),
    E0 = (31.2497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25868,0.0656455,-5.80652e-05,2.73511e-08,-5.50148e-12,3852.64,26.7202], Tmin=(100,'K'), Tmax=(1137.5,'K')), NASAPolynomial(coeffs=[10.0551,0.0347131,-1.72751e-05,3.44485e-09,-2.47352e-13,1851.45,-16.8508], Tmin=(1137.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.2497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C(C=C[O])=COO(22765)',
    structure = SMILES('C=C([CH]OO)C=C[O]'),
    E0 = (24.6649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0762264,0.0675569,-2.1069e-05,-3.83509e-08,2.4148e-11,3124.95,29.0628], Tmin=(100,'K'), Tmax=(942.661,'K')), NASAPolynomial(coeffs=[25.0332,0.00837739,-1.24397e-06,2.05828e-10,-2.13653e-14,-3656.07,-100.877], Tmin=(942.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.6649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C]1CC([O])C1O[O](22766)',
    structure = SMILES('[CH2][C]1CC([O])C1O[O]'),
    E0 = (485.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4338,0.0440355,5.75417e-06,-3.92134e-08,1.82817e-11,58552.7,32.3255], Tmin=(100,'K'), Tmax=(975.359,'K')), NASAPolynomial(coeffs=[13.5429,0.0238566,-8.55155e-06,1.55408e-09,-1.10759e-13,54788.2,-32.9803], Tmin=(975.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1C[CH]OC1O[O](22767)',
    structure = SMILES('[CH2][C]1C[CH]OC1O[O]'),
    E0 = (368.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5279,0.0428075,1.40422e-05,-5.79681e-08,2.98238e-11,44401.7,29.1132], Tmin=(100,'K'), Tmax=(865.884,'K')), NASAPolynomial(coeffs=[14.0961,0.0186256,-2.7538e-06,1.4846e-10,-1.58599e-15,40955.2,-37.0447], Tmin=(865.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOCs) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C[CH][O])C1OO1(22768)',
    structure = SMILES('[CH2][C](C[CH][O])C1OO1'),
    E0 = (502.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907714,0.0642149,-5.90472e-05,3.03523e-08,-6.35605e-12,60556.5,35.8183], Tmin=(100,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[11.8418,0.0261827,-9.43904e-06,1.59347e-09,-1.04012e-13,58041.7,-18.4601], Tmin=(1150,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(dioxirane) + radical(CCOJ) + radical(C2CJCOOH) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=[C]OO)C[CH][O](22769)',
    structure = SMILES('[CH2]C(=[C]OO)C[CH][O]'),
    E0 = (497.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.143729,0.0953472,-0.00015725,1.42003e-07,-4.98129e-11,59969.5,35.0999], Tmin=(100,'K'), Tmax=(845.06,'K')), NASAPolynomial(coeffs=[7.8599,0.03804,-1.86368e-05,3.54825e-09,-2.42847e-13,59407.5,3.56372], Tmin=(845.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'CC([CH][CH][O])=CO[O](22770)',
    structure = SMILES('C[C]([CH]O[O])C=C[O]'),
    E0 = (289.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.285127,0.0730391,-6.87997e-05,3.29521e-08,-6.21613e-12,34927.9,31.2528], Tmin=(100,'K'), Tmax=(1290.63,'K')), NASAPolynomial(coeffs=[17.2882,0.0203421,-7.55408e-06,1.31607e-09,-8.81213e-14,30539,-55.1155], Tmin=(1290.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(CCJ(C)CO) + radical(CCsJOOH)"""),
)

species(
    label = 'CC(=[C]O[O])C[CH][O](22771)',
    structure = SMILES('CC(=[C]O[O])C[CH][O]'),
    E0 = (498.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396392,0.0943992,-0.000171095,1.65505e-07,-5.99493e-11,60016.8,34.6707], Tmin=(100,'K'), Tmax=(872.431,'K')), NASAPolynomial(coeffs=[3.9861,0.0424203,-2.06541e-05,3.87868e-09,-2.61381e-13,60742.3,25.5897], Tmin=(872.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])CC[O](22772)',
    structure = SMILES('[CH2]C(=[C]O[O])CC[O]'),
    E0 = (469.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,180,2514.67,2515.21],'cm^-1')),
        HinderedRotor(inertia=(0.288251,'amu*angstrom^2'), symmetry=1, barrier=(6.62745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288198,'amu*angstrom^2'), symmetry=1, barrier=(6.62623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.33877,'amu*angstrom^2'), symmetry=1, barrier=(53.7728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288209,'amu*angstrom^2'), symmetry=1, barrier=(6.62648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651792,0.0833179,-0.000131827,1.20435e-07,-4.29826e-11,56549.1,34.4859], Tmin=(100,'K'), Tmax=(847.089,'K')), NASAPolynomial(coeffs=[5.68014,0.0395305,-1.87978e-05,3.54758e-09,-2.42083e-13,56416.3,15.3061], Tmin=(847.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C[CH]O(22773)',
    structure = SMILES('[CH2]C(=[C]O[O])C[CH]O'),
    E0 = (423.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,285.655,285.658],'cm^-1')),
        HinderedRotor(inertia=(0.00206577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119167,'amu*angstrom^2'), symmetry=1, barrier=(6.90033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119166,'amu*angstrom^2'), symmetry=1, barrier=(6.90034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119162,'amu*angstrom^2'), symmetry=1, barrier=(6.90033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119163,'amu*angstrom^2'), symmetry=1, barrier=(6.90033,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0506634,0.0944204,-0.000149235,1.25354e-07,-4.08861e-11,51111.6,35.3334], Tmin=(100,'K'), Tmax=(871.11,'K')), NASAPolynomial(coeffs=[10.9389,0.0300041,-1.34852e-05,2.46195e-09,-1.6407e-13,49761.8,-12.5536], Tmin=(871.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(C=CJO)"""),
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
    E0 = (378.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (378.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (558.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (588.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (794.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (608.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (792.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (993.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (861.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (870.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (386.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (386.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (441.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (752.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (609.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (503.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (452.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (442.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (488.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (500.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (531.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (396.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (378.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (480.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (642.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (805.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (671.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (837.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (490.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (554.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (653.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (537.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (489.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (432.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (471.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (506.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (681.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (973.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (386.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (385.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (386.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (461.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (433.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (441.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (403.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (503.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (452.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (502.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (639.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (580.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (633.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (513.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (456.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['vinoxy(1351)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C([CH]O[O])CC=O(21027)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]C(C[CH][O])O[O](21033)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', 'C=C([CH][O])C[CH][O](22742)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O2(2)', '[CH]C(=C)C[CH][O](15924)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][O](1548)', '[CH2]C(=C)[CH]O[O](21386)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]O[O](21387)', 'C=[C]C[CH][O](2773)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C=C([C]O[O])C[CH][O](22743)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C=C([CH]O[O])C[C][O](22744)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C(C[CH][O])C1OO1(22745)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C1CC([O])C1O[O](22746)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C(C=C[O])CO[O](10577)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH][CH][O])[CH]O[O](22747)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[O][CH]C[C]1CC1O[O](22748)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[O]O[CH][C]1CC([O])C1(22749)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[O]O[CH][C]1C[CH]OC1(22750)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[O][CH]C[C]1[CH]OOC1(22751)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.80789e+10,'s^-1'), n=0.280222, Ea=(64.2368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 60.6 to 64.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2]C1([CH]O[O])CC1[O](22752)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(110.296,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 109.0 to 110.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2]C1([CH]O[O])C[CH]O1(22753)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2]C1(C[CH][O])[CH]OO1(22754)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(153.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;doublebond_intra_2H;radadd_intra_O] + [R5;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'C=C([CH]O[O])C=C[O](21133)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O(4)', 'C=C(C=O)C[CH][O](22755)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(36.4151,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 32.2 to 36.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['vinoxy(1351)', '[CH2][C]=CO[O](16807)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][O](1556)', 'C=C=CO[O](16806)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00267131,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][O](1556)', '[CH2][C]=CO[O](16807)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'C=C([CH][CH][O])[CH]O[O](22756)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH]=C([CH]O[O])C[CH][O](20379)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C([CH]C[O])[CH]O[O](22757)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2][C]([CH]C=O)CO[O](10571)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00703183,'s^-1'), n=4.63833, Ea=(175.937,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeO;XH_out] + [R3H_SS_2Cd;C_rad_out_1H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(C[CH][O])CO[O](19320)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C([CH][CH]O)[CH]O[O](22758)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/Cd] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([CH]O[O])CC[O](20384)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C([CH]O[O])C[CH]O(20385)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['C=C([CH][CH][O])[CH]OO(22759)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(72901.1,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_2;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([CH]OO)C[CH][O](20381)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O][CH]CC[C]=CO[O](21034)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(T)(20)', '[O][CH]C[C]=CO[O](22760)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[O]OC=C1CC([O])C1(22761)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[O][CH]CC1=COOC1(19175)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2]C1=COOC([O])C1(22762)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSDSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['O(4)', '[O][CH]CC1=COC1(20375)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['O(4)', '[CH2]C1=COC([O])C1(22763)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_SSD;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['CC(C=C[O])=CO[O](22764)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2]C(C=C[O])=COO(22765)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2][C]1CC([O])C1O[O](22766)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2][C]1C[CH]OC1O[O](22767)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['[CH2][C](C[CH][O])C1OO1(22768)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(124.275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 121.7 to 124.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(=[C]OO)C[CH][O](22769)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['CC([CH][CH][O])=CO[O](22770)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C([CH]O[O])C[CH][O](20382)'],
    products = ['CC(=[C]O[O])C[CH][O](22771)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C(=[C]O[O])CC[O](22772)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C(=[C]O[O])C[CH]O(22773)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4710',
    isomers = [
        'C=C([CH]O[O])C[CH][O](20382)',
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
    label = 'PDepNetwork #4710',
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


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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.644475,0.0716765,-6.98211e-05,3.5401e-08,-7.22043e-12,86995.4,32.3759], Tmin=(100,'K'), Tmax=(1178.73,'K')), NASAPolynomial(coeffs=[14.3121,0.0252958,-1.07995e-05,2.01971e-09,-1.40535e-13,83773.3,-35.81], Tmin=(1178.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]C(=C)C(=C)[CH]O[O](20338)',
    structure = SMILES('[CH]C(=C)C(=C)[CH]O[O]'),
    E0 = (577.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123104,0.0750846,-5.47386e-05,1.21598e-08,2.17045e-12,69556.9,30.0934], Tmin=(100,'K'), Tmax=(1003.87,'K')), NASAPolynomial(coeffs=[17.7075,0.0248468,-9.30081e-06,1.66062e-09,-1.14836e-13,65027.2,-59.7856], Tmin=(1003.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC([C]=C)O[O](21285)',
    structure = SMILES('[CH]=[C]CC([C]=C)O[O]'),
    E0 = (837.183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,309.244,309.502],'cm^-1')),
        HinderedRotor(inertia=(0.15832,'amu*angstrom^2'), symmetry=1, barrier=(10.7595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15865,'amu*angstrom^2'), symmetry=1, barrier=(10.7609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158659,'amu*angstrom^2'), symmetry=1, barrier=(10.7634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158627,'amu*angstrom^2'), symmetry=1, barrier=(10.7623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3940.4,'J/mol'), sigma=(6.57992,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.48 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556652,0.080901,-0.000105444,7.68688e-08,-2.28706e-11,100809,33.5511], Tmin=(100,'K'), Tmax=(816.381,'K')), NASAPolynomial(coeffs=[10.6219,0.0315779,-1.48073e-05,2.84365e-09,-1.98929e-13,99165.9,-12.9648], Tmin=(816.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(837.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC(=C)[CH][O](22626)',
    structure = SMILES('[CH]=[C]CC([CH2])=C[O]'),
    E0 = (618.802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,274.114,274.206],'cm^-1')),
        HinderedRotor(inertia=(0.42717,'amu*angstrom^2'), symmetry=1, barrier=(22.8358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426678,'amu*angstrom^2'), symmetry=1, barrier=(22.8308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.424746,'amu*angstrom^2'), symmetry=1, barrier=(22.8339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676507,0.0622978,-4.00584e-05,-2.07606e-09,7.90482e-12,74554.1,27.9358], Tmin=(100,'K'), Tmax=(957.817,'K')), NASAPolynomial(coeffs=[18.2745,0.0141735,-4.42034e-06,7.70382e-10,-5.54581e-14,70019.3,-62.2804], Tmin=(957.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC(=[CH])[CH2](16901)',
    structure = SMILES('[CH]C(=C)C[C]=[CH]'),
    E0 = (905.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,351.264,352.028,352.055,352.809],'cm^-1')),
        HinderedRotor(inertia=(0.584914,'amu*angstrom^2'), symmetry=1, barrier=(51.0487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585247,'amu*angstrom^2'), symmetry=1, barrier=(51.0279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58312,'amu*angstrom^2'), symmetry=1, barrier=(51.0421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.68,'J/mol'), sigma=(5.67268,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.00 K, Pc=40.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46847,0.0539598,-4.00949e-05,1.60924e-08,-2.72863e-12,108977,24.6512], Tmin=(100,'K'), Tmax=(1348.6,'K')), NASAPolynomial(coeffs=[10.0074,0.028633,-1.19247e-05,2.16667e-09,-1.47114e-13,106674,-19.0978], Tmin=(1348.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=[CH](9646)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([554.803,1738.79,3454.47],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97621,0.00212915,-8.08978e-08,-3.83305e-10,9.76908e-14,101881,6.00119], Tmin=(100,'K'), Tmax=(1982.31,'K')), NASAPolynomial(coeffs=[5.05695,0.00131032,-4.91873e-07,1.01502e-10,-7.16167e-15,101185,-0.627232], Tmin=(1982.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
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
    label = '[CH]=[C]C[C]=C(15982)',
    structure = SMILES('[CH]=[C]C[C]=C'),
    E0 = (811.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,250.816],'cm^-1')),
        HinderedRotor(inertia=(0.246985,'amu*angstrom^2'), symmetry=1, barrier=(11.0542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248373,'amu*angstrom^2'), symmetry=1, barrier=(11.0546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3393,0.0378852,-3.20383e-05,1.55074e-08,-3.23394e-12,97663.1,20.1382], Tmin=(100,'K'), Tmax=(1106.63,'K')), NASAPolynomial(coeffs=[7.00928,0.0210051,-9.15799e-06,1.72363e-09,-1.2002e-13,96629.5,-2.86492], Tmin=(1106.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(811.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC(=C)[C]O[O](22627)',
    structure = SMILES('[CH]=[C]CC([CH2])=[C]O[O]'),
    E0 = (993.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,345.109],'cm^-1')),
        HinderedRotor(inertia=(0.0962577,'amu*angstrom^2'), symmetry=1, barrier=(8.13515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962577,'amu*angstrom^2'), symmetry=1, barrier=(8.13517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962565,'amu*angstrom^2'), symmetry=1, barrier=(8.13514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962559,'amu*angstrom^2'), symmetry=1, barrier=(8.13516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.527548,0.0818049,-0.000120922,9.82231e-08,-3.17186e-11,119620,34.9925], Tmin=(100,'K'), Tmax=(834.404,'K')), NASAPolynomial(coeffs=[10.3124,0.0282924,-1.28485e-05,2.38716e-09,-1.62044e-13,118217,-9.0644], Tmin=(834.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(993.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]CC(=C)[CH]O[O](22628)',
    structure = SMILES('[C]=[C]CC(=C)[CH]O[O]'),
    E0 = (1033.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,409.786,409.792,409.837],'cm^-1')),
        HinderedRotor(inertia=(0.0920541,'amu*angstrom^2'), symmetry=1, barrier=(10.9713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0920501,'amu*angstrom^2'), symmetry=1, barrier=(10.971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00100373,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347226,'amu*angstrom^2'), symmetry=1, barrier=(41.3938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852681,0.0723791,-8.22858e-05,5.03673e-08,-1.25555e-11,124389,31.6136], Tmin=(100,'K'), Tmax=(967.089,'K')), NASAPolynomial(coeffs=[11.6822,0.0275866,-1.28106e-05,2.47415e-09,-1.74721e-13,122294,-20.2703], Tmin=(967.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1033.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC(=C)C1OO1(22629)',
    structure = SMILES('[CH]=[C]CC(=C)C1OO1'),
    E0 = (568.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.608903,0.0648304,-4.22336e-05,3.189e-10,7.09623e-12,68540,30.0399], Tmin=(100,'K'), Tmax=(944.921,'K')), NASAPolynomial(coeffs=[17.2143,0.0181516,-5.62054e-06,9.35136e-10,-6.41465e-14,64347.6,-54.7092], Tmin=(944.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC(=C)C1O[O](22630)',
    structure = SMILES('[CH]=C1CC(=C)C1O[O]'),
    E0 = (488.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11738,0.0486016,2.32781e-06,-4.07064e-08,1.94824e-11,58919.1,27.6786], Tmin=(100,'K'), Tmax=(989.803,'K')), NASAPolynomial(coeffs=[16.5997,0.0206773,-7.85421e-06,1.51209e-09,-1.12173e-13,54157.2,-55.4284], Tmin=(989.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC(=C)CO[O](20303)',
    structure = SMILES('C#CC=C([CH2])CO[O]'),
    E0 = (387.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69524,0.0758728,-8.93989e-05,5.97799e-08,-1.63808e-11,46711.7,26.2492], Tmin=(100,'K'), Tmax=(883.404,'K')), NASAPolynomial(coeffs=[10.5692,0.0311665,-1.34929e-05,2.50007e-09,-1.71721e-13,44967.1,-20.1634], Tmin=(883.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(ROOJ) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=[C][CH]C([CH2])[CH]O[O](22631)',
    structure = SMILES('[CH][C]=CC([CH2])[CH]O[O]'),
    E0 = (981.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.238773,0.089655,-0.000126265,1.04626e-07,-3.46221e-11,118185,35.1224], Tmin=(100,'K'), Tmax=(866.502,'K')), NASAPolynomial(coeffs=[8.04331,0.0402474,-1.75739e-05,3.1809e-09,-2.1215e-13,117335,1.48679], Tmin=(866.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(981.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[C]1CC1O[O](22632)',
    structure = SMILES('[CH]=[C]C[C]1CC1O[O]'),
    E0 = (816.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13546,0.0541984,-2.9447e-05,7.41971e-10,3.13945e-12,98296.1,33.7263], Tmin=(100,'K'), Tmax=(1084.53,'K')), NASAPolynomial(coeffs=[13.2831,0.0247662,-9.99973e-06,1.85638e-09,-1.29987e-13,94757.2,-30.0327], Tmin=(1084.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(816.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(C2CJCOOH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[C]([CH]O[O])C1(22633)',
    structure = SMILES('[CH]=C1C[C]([CH]O[O])C1'),
    E0 = (773.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5416,0.0520407,-3.12634e-05,8.96912e-09,-1.03835e-12,93115.7,30.6797], Tmin=(100,'K'), Tmax=(1917.66,'K')), NASAPolynomial(coeffs=[13.7481,0.0265789,-1.13466e-05,2.04495e-09,-1.35639e-13,88434.3,-36.1571], Tmin=(1917.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(773.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOOH) + radical(Cds_P)"""),
)

species(
    label = '[O]O[CH][C]1C[C]=CC1(22634)',
    structure = SMILES('[O]O[CH][C]1C[C]=CC1'),
    E0 = (689.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60208,0.0441506,-8.73196e-06,-1.34757e-08,6.36418e-12,83064.1,30.9357], Tmin=(100,'K'), Tmax=(1125.01,'K')), NASAPolynomial(coeffs=[10.417,0.0299427,-1.26331e-05,2.37367e-09,-1.66217e-13,79996.5,-17.4485], Tmin=(1125.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOOH) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]=[C]C[C]1[CH]OOC1(22635)',
    structure = SMILES('[CH]=[C]C[C]1[CH]OOC1'),
    E0 = (785.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46372,0.0451412,5.25402e-06,-4.39881e-08,2.27343e-11,94547.4,31.0022], Tmin=(100,'K'), Tmax=(893.583,'K')), NASAPolynomial(coeffs=[13.0966,0.0223429,-5.61695e-06,7.84517e-10,-4.8872e-14,91299.6,-30.3505], Tmin=(893.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(C2CJCOOH) + radical(CCsJOOC) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC1([CH2])[CH]O[O](22636)',
    structure = SMILES('[CH]=C1CC1([CH2])[CH]O[O]'),
    E0 = (826.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576635,0.0745654,-7.49275e-05,3.91594e-08,-8.24611e-12,99570.4,28.444], Tmin=(100,'K'), Tmax=(1141.47,'K')), NASAPolynomial(coeffs=[14.3107,0.0264372,-1.16819e-05,2.22094e-09,-1.55943e-13,96435,-39.6322], Tmin=(1141.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(826.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(CCsJOOH) + radical(Neopentyl) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC1([CH2])[CH]OO1(22638)',
    structure = SMILES('[CH]=[C]CC1([CH2])[CH]OO1'),
    E0 = (872.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439376,0.0833005,-0.000108699,7.97761e-08,-2.3877e-11,105082,27.8401], Tmin=(100,'K'), Tmax=(812.193,'K')), NASAPolynomial(coeffs=[10.7465,0.0325423,-1.49627e-05,2.84046e-09,-1.97198e-13,103407,-19.7423], Tmin=(812.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(872.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(CJCOOH) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC(=C)C=O(22640)',
    structure = SMILES('[CH]=[C]CC(=C)C=O'),
    E0 = (442.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,310.177],'cm^-1')),
        HinderedRotor(inertia=(0.182809,'amu*angstrom^2'), symmetry=1, barrier=(12.4873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18295,'amu*angstrom^2'), symmetry=1, barrier=(12.4888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182882,'amu*angstrom^2'), symmetry=1, barrier=(12.4879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68419,0.0569697,-5.07711e-05,2.53406e-08,-5.58419e-12,53344.7,24.6538], Tmin=(100,'K'), Tmax=(1027.48,'K')), NASAPolynomial(coeffs=[7.63505,0.0338027,-1.69496e-05,3.39578e-09,-2.44662e-13,52121.9,-4.21697], Tmin=(1027.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC(=[CH])[CH]O[O](20305)',
    structure = SMILES('[CH]=[C]CC(=[CH])[CH]O[O]'),
    E0 = (969.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,1685,370,3115,3125,620,680,785,800,1600,1700,759.648],'cm^-1')),
        HinderedRotor(inertia=(0.549553,'amu*angstrom^2'), symmetry=1, barrier=(12.6353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206127,'amu*angstrom^2'), symmetry=1, barrier=(12.6804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549378,'amu*angstrom^2'), symmetry=1, barrier=(12.6313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097436,'amu*angstrom^2'), symmetry=1, barrier=(40.2899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643494,0.0742918,-8.31326e-05,4.8306e-08,-1.12062e-11,116712,32.7867], Tmin=(100,'K'), Tmax=(1046.67,'K')), NASAPolynomial(coeffs=[13.9271,0.0235261,-1.0379e-05,1.966e-09,-1.37708e-13,113931,-31.905], Tmin=(1046.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(969.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]C(=C)[CH]O[O](22642)',
    structure = SMILES('[CH]=CC=C([CH2])[CH]O[O]'),
    E0 = (571.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521768,0.0654208,-3.88406e-05,-3.9578e-09,8.38535e-12,68890.2,30.8758], Tmin=(100,'K'), Tmax=(963.466,'K')), NASAPolynomial(coeffs=[17.8649,0.0189882,-6.36045e-06,1.11406e-09,-7.8424e-14,64361.5,-58.3083], Tmin=(963.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(=C)CO[O](20241)',
    structure = SMILES('[CH]=[C]C=C([CH2])CO[O]'),
    E0 = (653.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.487694,0.0798819,-9.61005e-05,6.28509e-08,-1.59225e-11,78705.4,30.0088], Tmin=(100,'K'), Tmax=(775.404,'K')), NASAPolynomial(coeffs=[11.7158,0.0290811,-1.16018e-05,2.04425e-09,-1.35828e-13,76750.1,-22.6843], Tmin=(775.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(=[CH])CO[O](19501)',
    structure = SMILES('[CH]=[C]CC(=[CH])CO[O]'),
    E0 = (852.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,1685,370,3115,3125,620,680,785,800,1600,1700,229.602],'cm^-1')),
        HinderedRotor(inertia=(0.234971,'amu*angstrom^2'), symmetry=1, barrier=(8.79349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235104,'amu*angstrom^2'), symmetry=1, barrier=(8.79558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.23498,'amu*angstrom^2'), symmetry=1, barrier=(8.79398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235031,'amu*angstrom^2'), symmetry=1, barrier=(8.79457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3938.88,'J/mol'), sigma=(6.5717,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=615.24 K, Pc=31.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.361329,0.0863095,-0.000126995,1.04932e-07,-3.46822e-11,102609,33.1973], Tmin=(100,'K'), Tmax=(824.125,'K')), NASAPolynomial(coeffs=[9.77907,0.0327227,-1.51248e-05,2.83947e-09,-1.94131e-13,101325,-8.79339], Tmin=(824.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(852.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

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
    label = '[CH]=CCC(=[CH])[CH]O[O](20308)',
    structure = SMILES('[CH]=CCC(=[CH])[CH]O[O]'),
    E0 = (731.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3115,3125,620,680,785,800,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(0.0763635,'amu*angstrom^2'), symmetry=1, barrier=(1.75575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771838,'amu*angstrom^2'), symmetry=1, barrier=(17.7461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771846,'amu*angstrom^2'), symmetry=1, barrier=(17.7462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883439,'amu*angstrom^2'), symmetry=1, barrier=(31.3646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400592,0.0729746,-6.98807e-05,3.39352e-08,-6.52091e-12,88120.6,33.1813], Tmin=(100,'K'), Tmax=(1261.49,'K')), NASAPolynomial(coeffs=[16.6403,0.0214813,-8.65228e-06,1.57781e-09,-1.08439e-13,84023.4,-48.9385], Tmin=(1261.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(731.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[CH]=[C][CH]C(=C)[CH]OO(22643)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH]OO'),
    E0 = (618.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.143002,0.0759871,-6.47094e-05,2.14744e-08,-3.88103e-13,74553.4,31.073], Tmin=(100,'K'), Tmax=(977.763,'K')), NASAPolynomial(coeffs=[18.2412,0.0201728,-7.0433e-06,1.21952e-09,-8.34772e-14,70143.1,-60.2885], Tmin=(977.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(=[CH])[CH]OO(20306)',
    structure = SMILES('[CH]=[C]CC(=[CH])[CH]OO'),
    E0 = (817.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,1685,370,3115,3125,620,680,785,800,1600,1700,1395.73],'cm^-1')),
        HinderedRotor(inertia=(0.644249,'amu*angstrom^2'), symmetry=1, barrier=(21.6161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64432,'amu*angstrom^2'), symmetry=1, barrier=(21.6161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644301,'amu*angstrom^2'), symmetry=1, barrier=(21.6161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86708,'amu*angstrom^2'), symmetry=1, barrier=(96.1716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103404,'amu*angstrom^2'), symmetry=1, barrier=(21.6161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.325379,0.0785533,-8.08465e-05,4.21182e-08,-8.72396e-12,98443.9,33.1666], Tmin=(100,'K'), Tmax=(1168.07,'K')), NASAPolynomial(coeffs=[16.3005,0.0238476,-1.05956e-05,2.02332e-09,-1.42592e-13,94711.9,-46.386], Tmin=(1168.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(817.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC[C]=CO[O](21286)',
    structure = SMILES('[CH]=[C]CC[C]=CO[O]'),
    E0 = (854.187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.387146,'amu*angstrom^2'), symmetry=1, barrier=(8.90125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387124,'amu*angstrom^2'), symmetry=1, barrier=(8.90075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.38714,'amu*angstrom^2'), symmetry=1, barrier=(8.9011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387358,'amu*angstrom^2'), symmetry=1, barrier=(8.90612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3996.75,'J/mol'), sigma=(6.61607,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=624.28 K, Pc=31.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.466307,0.082176,-0.00011149,8.39885e-08,-2.5519e-11,102858,33.742], Tmin=(100,'K'), Tmax=(804.059,'K')), NASAPolynomial(coeffs=[11.0634,0.0294593,-1.31476e-05,2.45189e-09,-1.68014e-13,101154,-15.0719], Tmin=(804.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(854.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]C[C]=CO[O](22448)',
    structure = SMILES('[CH]=[C]C[C]=CO[O]'),
    E0 = (879.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.463323,'amu*angstrom^2'), symmetry=1, barrier=(10.6527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.463752,'amu*angstrom^2'), symmetry=1, barrier=(10.6626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.463145,'amu*angstrom^2'), symmetry=1, barrier=(10.6486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09168,0.0675752,-9.83476e-05,7.69228e-08,-2.3847e-11,105849,29.2755], Tmin=(100,'K'), Tmax=(833.866,'K')), NASAPolynomial(coeffs=[10.2962,0.0206986,-9.1256e-06,1.67487e-09,-1.12998e-13,104409,-12.8913], Tmin=(833.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(879.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC(=CO[O])C1(22644)',
    structure = SMILES('[CH]=C1CC(=CO[O])C1'),
    E0 = (510.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15709,0.045875,1.36905e-05,-5.7501e-08,2.71698e-11,61569.3,28.2908], Tmin=(100,'K'), Tmax=(952.536,'K')), NASAPolynomial(coeffs=[17.9692,0.016465,-4.85897e-06,8.78029e-10,-6.68846e-14,56497.8,-61.8092], Tmin=(952.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC1=COOC1(19220)',
    structure = SMILES('[CH]=[C]CC1=COOC1'),
    E0 = (534.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925413,0.0567049,-2.61427e-05,-8.69893e-09,7.72605e-12,64364.8,29.8319], Tmin=(100,'K'), Tmax=(1021.56,'K')), NASAPolynomial(coeffs=[15.5342,0.0221651,-8.70205e-06,1.63467e-09,-1.17033e-13,60197.6,-46.7466], Tmin=(1021.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolene) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC([CH2])=COO1(22645)',
    structure = SMILES('[CH]=C1CC([CH2])=COO1'),
    E0 = (422.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34385,0.037029,4.80892e-05,-8.87028e-08,3.54704e-11,50891.7,23.4642], Tmin=(100,'K'), Tmax=(1002.05,'K')), NASAPolynomial(coeffs=[17.7379,0.0238006,-1.02692e-05,2.12373e-09,-1.63325e-13,44984.8,-68.7414], Tmin=(1002.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC1=COC1(20301)',
    structure = SMILES('[CH]=[C]CC1=COC1'),
    E0 = (533.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10635,0.0446545,1.85422e-05,-6.90707e-08,3.33169e-11,64266.3,23.7229], Tmin=(100,'K'), Tmax=(934.21,'K')), NASAPolynomial(coeffs=[20.8792,0.00833302,-7.55243e-07,8.85956e-11,-1.27668e-14,58462.4,-81.6142], Tmin=(934.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC([CH2])=CO1(22646)',
    structure = SMILES('[CH]=C1CC(=C)[CH]O1'),
    E0 = (304.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03553,0.0191901,8.40268e-05,-1.29729e-07,5.28853e-11,36772.2,22.5084], Tmin=(100,'K'), Tmax=(940.775,'K')), NASAPolynomial(coeffs=[18.0531,0.0121789,-2.20206e-06,4.02312e-10,-3.86471e-14,31054.9,-68.1578], Tmin=(940.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJ(O)C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC(C)=CO[O](22647)',
    structure = SMILES('C#CC=C(C)[CH]O[O]'),
    E0 = (383.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.333577,0.0717868,-6.51347e-05,2.96585e-08,-5.32853e-12,46239.2,28.4915], Tmin=(100,'K'), Tmax=(1348.05,'K')), NASAPolynomial(coeffs=[17.5999,0.0205526,-8.1248e-06,1.46438e-09,-9.97873e-14,41584.1,-59.9653], Tmin=(1348.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH]=C=CC([CH2])=COO(22648)',
    structure = SMILES('C#CC=C([CH2])[CH]OO'),
    E0 = (352.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.239603,0.073215,-6.18784e-05,2.24967e-08,-2.04117e-12,42564.5,27.7149], Tmin=(100,'K'), Tmax=(1059.67,'K')), NASAPolynomial(coeffs=[17.906,0.0208957,-8.15556e-06,1.49257e-09,-1.04285e-13,38013.7,-62.3451], Tmin=(1059.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CCJO) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=C1C[C]([CH2])C1O[O](22649)',
    structure = SMILES('[CH]=C1C[C]([CH2])C1O[O]'),
    E0 = (776.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3761,0.0483482,-1.38677e-05,-1.37148e-08,7.88579e-12,93552.2,33.4672], Tmin=(100,'K'), Tmax=(1041.96,'K')), NASAPolynomial(coeffs=[12.0181,0.0267838,-1.05926e-05,1.95681e-09,-1.37254e-13,90287.4,-23.3363], Tmin=(1041.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(776.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]1C[C]=CC1O[O](22650)',
    structure = SMILES('[CH2][C]1C[C]=CC1O[O]'),
    E0 = (693.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70845,0.0372907,1.97265e-05,-5.06845e-08,2.15955e-11,83488.3,32.0497], Tmin=(100,'K'), Tmax=(980.5,'K')), NASAPolynomial(coeffs=[12.5754,0.0245108,-8.9924e-06,1.66211e-09,-1.19578e-13,79840.6,-27.8974], Tmin=(980.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(693.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH]=[C]C[C]([CH2])C1OO1(22651)',
    structure = SMILES('[CH]=[C]C[C]([CH2])C1OO1'),
    E0 = (841.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07038,0.0538941,-1.70116e-05,-2.35156e-08,1.55509e-11,101380,35.718], Tmin=(100,'K'), Tmax=(917.186,'K')), NASAPolynomial(coeffs=[15.4527,0.0188895,-5.09651e-06,7.74397e-10,-5.13001e-14,97575.5,-38.781], Tmin=(917.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(841.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(C2CJCOOH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC([CH2])=[C]OO(22652)',
    structure = SMILES('[CH]=[C]CC([CH2])=[C]OO'),
    E0 = (841.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,298.231],'cm^-1')),
        HinderedRotor(inertia=(0.152337,'amu*angstrom^2'), symmetry=1, barrier=(9.61475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152336,'amu*angstrom^2'), symmetry=1, barrier=(9.61475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00189538,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311169,'amu*angstrom^2'), symmetry=1, barrier=(19.6393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871061,'amu*angstrom^2'), symmetry=1, barrier=(54.9767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437248,0.0831711,-0.00010756,7.63619e-08,-2.1978e-11,101342,34.5694], Tmin=(100,'K'), Tmax=(844.966,'K')), NASAPolynomial(coeffs=[11.5439,0.0305934,-1.42233e-05,2.72106e-09,-1.90059e-13,99465.5,-17.1428], Tmin=(844.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(841.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(C)=CO[O](22653)',
    structure = SMILES('[CH]=C=C[C](C)[CH]O[O]'),
    E0 = (651.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666102,0.0763,-9.05373e-05,6.14582e-08,-1.70849e-11,78489.2,30.1063], Tmin=(100,'K'), Tmax=(872.114,'K')), NASAPolynomial(coeffs=[10.4526,0.0314126,-1.33312e-05,2.43854e-09,-1.66009e-13,76782.2,-15.7687], Tmin=(872.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCJ(C)CO) + radical(CCsJOOH) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]CC(C)=[C]O[O](22654)',
    structure = SMILES('[CH]=[C]CC(C)=[C]O[O]'),
    E0 = (842.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.515554,0.084474,-0.000130266,1.1294e-07,-3.84879e-11,101397,34.7534], Tmin=(100,'K'), Tmax=(850.816,'K')), NASAPolynomial(coeffs=[8.06916,0.0342583,-1.58125e-05,2.94758e-09,-1.99799e-13,100644,2.65944], Tmin=(850.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(842.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC([CH2])=[C]O[O](22655)',
    structure = SMILES('[CH]=CCC([CH2])=[C]O[O]'),
    E0 = (755.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,266.405],'cm^-1')),
        HinderedRotor(inertia=(0.00235639,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230916,'amu*angstrom^2'), symmetry=1, barrier=(11.6388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18868,'amu*angstrom^2'), symmetry=1, barrier=(9.52172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189338,'amu*angstrom^2'), symmetry=1, barrier=(9.50556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726407,0.0751073,-8.80794e-05,5.73527e-08,-1.52251e-11,91009.9,33.8143], Tmin=(100,'K'), Tmax=(911.953,'K')), NASAPolynomial(coeffs=[11.14,0.0294308,-1.29488e-05,2.4291e-09,-1.68364e-13,89110.6,-15.4652], Tmin=(911.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (722.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (967.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (931.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1138.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (952.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1208.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1337.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1205.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1245.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (730.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (730.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (785.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1004.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (953.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (848.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (754.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (785.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (827.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (789.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (874.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (728.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (722.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (802.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (889.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1051.8,'kJ/mol'),
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
    E0 = (1181.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (839.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (898.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (997.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (766.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (775.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (913.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (815.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (850.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1024.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (1317.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (730.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (729.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (730.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (805.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (777.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (785.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (747.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (848.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (754.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (841.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (983.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (924.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (977.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (800.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (928.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['C3H3(5450)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]C(=C)C(=C)[CH]O[O](20338)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC([C]=C)O[O](21285)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]CC(=C)[CH][O](22626)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O2(2)', '[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=[CH](9646)', '[CH2]C(=C)[CH]O[O](21386)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]O[O](21387)', '[CH]=[C]C[C]=C(15982)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]=[C]CC(=C)[C]O[O](22627)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[C]=[C]CC(=C)[CH]O[O](22628)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]CC(=C)C1OO1(22629)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C1CC(=C)C1O[O](22630)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C=CC(=C)CO[O](20303)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C][CH]C([CH2])[CH]O[O](22631)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]C[C]1CC1O[O](22632)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C1C[C]([CH]O[O])C1(22633)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[O]O[CH][C]1C[C]=CC1(22634)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]C[C]1[CH]OOC1(22635)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.80789e+10,'s^-1'), n=0.280222, Ea=(62.9673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 59.4 to 63.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C1CC1([CH2])[CH]O[O](22636)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(105.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH2]C1([CH]O[O])C=[C]C1(22637)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(67.4564,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 66.6 to 67.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]CC1([CH2])[CH]OO1(22638)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(151.972,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;doublebond_intra_2H;radadd_intra_O] + [R5;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH]=C=CC(=C)[CH]O[O](22639)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH]=[C]CC(=C)C=O(22640)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(36.4151,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 32.2 to 36.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C3H3(5450)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0130473,'m^3/(mol*s)'), n=2.32138, Ea=(27.0133,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C][CH2](16918)', 'C=C=CO[O](16806)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00267131,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C][CH2](16918)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH]=[C][CH]C(=C)[CH]O[O](22641)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH]=[C]CC(=[CH])[CH]O[O](20305)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C[CH]C(=C)[CH]O[O](22642)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.18474e+10,'s^-1'), n=0.962222, Ea=(117.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C][CH]C(=C)CO[O](20241)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00703183,'s^-1'), n=4.63833, Ea=(175.937,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeO;XH_out] + [R3H_SS_2Cd;C_rad_out_1H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]CC(=[CH])CO[O](19501)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CCC(=[CH])[CH]O[O](20308)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([CH]O[O])C[C]=C(20080)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C][CH]C(=C)[CH]OO(22643)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(72901.1,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_2;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C]CC(=[CH])[CH]OO(20306)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC[C]=CO[O](21286)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(T)(20)', '[CH]=[C]C[C]=CO[O](22448)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C1CC(=CO[O])C1(22644)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]CC1=COOC1(19220)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C1CC([CH2])=COO1(22645)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSDSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['O(4)', '[CH]=[C]CC1=COC1(20301)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['O(4)', '[CH]=C1CC([CH2])=CO1(22646)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_SSD;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C=CC(C)=CO[O](22647)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C=CC([CH2])=COO(22648)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=C1C[C]([CH2])C1O[O](22649)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH2][C]1C[C]=CC1O[O](22650)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]C[C]([CH2])C1OO1(22651)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(119.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 117.1 to 119.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=[C]CC([CH2])=[C]OO(22652)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C][CH]C(C)=CO[O](22653)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    products = ['[CH]=[C]CC(C)=[C]O[O](22654)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=CCC([CH2])=[C]O[O](22655)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=[C]CC(=C)[C]O[O](22313)'],
    products = ['[CH]=[C]CC(=C)[CH]O[O](20307)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #5004',
    isomers = [
        '[CH]=[C]CC(=C)[CH]O[O](20307)',
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
    label = 'PDepNetwork #5004',
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


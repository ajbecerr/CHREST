species(
    label = '[CH2][C]([O])[CH]OO[C]=C(25548)',
    structure = SMILES('[CH2][C]([O])[CH]OO[C]=C'),
    E0 = (768.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0565599,0.10018,-0.000172155,1.5472e-07,-5.35202e-11,92593.6,35.799], Tmin=(100,'K'), Tmax=(850.351,'K')), NASAPolynomial(coeffs=[9.4802,0.0337623,-1.69688e-05,3.24337e-09,-2.21663e-13,91751.1,-4.08157], Tmin=(850.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(768.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOOC) + radical(CJCO) + radical(C=CJO)"""),
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
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4007.73,'J/mol'), sigma=(6.4029,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=626.00 K, Pc=34.64 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])=COO[C]=C(25530)',
    structure = SMILES('C=[C]OO[CH]C(=C)[O]'),
    E0 = (310.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.468817,0.0690217,-6.49693e-05,3.02761e-08,-5.52827e-12,37427.5,34.4151], Tmin=(100,'K'), Tmax=(1332.74,'K')), NASAPolynomial(coeffs=[17.6348,0.0175009,-6.98264e-06,1.2699e-09,-8.71854e-14,32851.9,-53.3314], Tmin=(1332.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]O[O](1354)',
    structure = SMILES('C=[C]O[O]'),
    E0 = (337.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.00100972,'amu*angstrom^2'), symmetry=1, barrier=(11.4644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.15943,0.0239683,-4.52732e-05,4.81856e-08,-1.85891e-11,40635.7,16.5558], Tmin=(100,'K'), Tmax=(873.282,'K')), NASAPolynomial(coeffs=[2.07161,0.0153864,-7.23324e-06,1.359e-09,-9.18757e-14,41343,24.6179], Tmin=(873.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH][C]([CH2])[O](13252)',
    structure = SMILES('[CH][C]([CH2])[O]'),
    E0 = (775.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,418.116,418.118,1773.98,1774.02],'cm^-1')),
        HinderedRotor(inertia=(0.000964274,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.067747,'amu*angstrom^2'), symmetry=1, barrier=(8.40473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31074,0.0403608,-6.035e-05,4.62802e-08,-1.34634e-11,93383.3,16.9577], Tmin=(100,'K'), Tmax=(685.403,'K')), NASAPolynomial(coeffs=[8.03541,0.0109907,-4.91293e-06,9.56161e-10,-6.74028e-14,92503.7,-9.19014], Tmin=(685.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(775.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][C]([O])[CH]O[O](22833)',
    structure = SMILES('[CH2][C]([O])[CH]O[O]'),
    E0 = (577.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,2362.58],'cm^-1')),
        HinderedRotor(inertia=(0.284489,'amu*angstrom^2'), symmetry=1, barrier=(6.54095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285315,'amu*angstrom^2'), symmetry=1, barrier=(6.55994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.59577,'amu*angstrom^2'), symmetry=1, barrier=(59.6819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05836,0.0757261,-0.00014534,1.35191e-07,-4.65652e-11,69594.6,24.5151], Tmin=(100,'K'), Tmax=(893.702,'K')), NASAPolynomial(coeffs=[7.28479,0.0216789,-1.06872e-05,1.96768e-09,-1.29125e-13,69527.1,1.02502], Tmin=(893.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH) + radical(CCsJOOH) + radical(CJCO)"""),
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
    label = '[CH2][C]([O])[C]OO[C]=C(26854)',
    structure = SMILES('[CH2][C]([O])[C]OO[C]=C'),
    E0 = (1043.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,360,370,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287585,0.104398,-0.000183522,1.6108e-07,-5.4266e-11,125623,35.0654], Tmin=(100,'K'), Tmax=(847.618,'K')), NASAPolynomial(coeffs=[12.4253,0.0264678,-1.38708e-05,2.68059e-09,-1.83764e-13,124112,-20.3644], Tmin=(847.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1043.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH][C]([O])[CH]OO[C]=C(26855)',
    structure = SMILES('[CH][C]([O])[CH]OO[C]=C'),
    E0 = (1005.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,360,370,350,2950,3100,1380,975,1025,1650,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0423117,0.0980582,-0.00017196,1.55279e-07,-5.37496e-11,121050,35.8144], Tmin=(100,'K'), Tmax=(851.573,'K')), NASAPolynomial(coeffs=[9.55883,0.0314706,-1.61175e-05,3.09356e-09,-2.1157e-13,120222,-3.90876], Tmin=(851.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCsJOOC) + radical(C2CsJOH) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]([O])C1OOC1=C(26856)',
    structure = SMILES('[CH2][C]([O])C1OOC1=C'),
    E0 = (470.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659426,0.0623889,-4.16271e-05,5.41363e-09,2.91352e-12,56763.1,29.6773], Tmin=(100,'K'), Tmax=(1077.51,'K')), NASAPolynomial(coeffs=[17.6561,0.0188602,-8.27092e-06,1.62955e-09,-1.18867e-13,51964.3,-58.8619], Tmin=(1077.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH]=[C]OOC=C([CH2])[O](25806)',
    structure = SMILES('[CH]=[C]OO[CH]C(=C)[O]'),
    E0 = (557.169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,555.289,555.293],'cm^-1')),
        HinderedRotor(inertia=(0.00965225,'amu*angstrom^2'), symmetry=1, barrier=(2.11208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100584,'amu*angstrom^2'), symmetry=1, barrier=(22.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957208,'amu*angstrom^2'), symmetry=1, barrier=(22.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7372,'amu*angstrom^2'), symmetry=1, barrier=(39.9416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.588317,0.0703458,-7.433e-05,3.87769e-08,-7.90956e-12,67138.7,34.3845], Tmin=(100,'K'), Tmax=(1198.03,'K')), NASAPolynomial(coeffs=[16.6924,0.0165771,-7.00862e-06,1.3146e-09,-9.2085e-14,63280.1,-46.2181], Tmin=(1198.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH2][C]([O])[CH][O](26240)',
    structure = SMILES('[CH2][C]([O])[CH][O]'),
    E0 = (571.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,180,2035.01],'cm^-1')),
        HinderedRotor(inertia=(0.0866393,'amu*angstrom^2'), symmetry=1, barrier=(1.99201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086712,'amu*angstrom^2'), symmetry=1, barrier=(1.99368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60835,0.065059,-0.00013406,1.3209e-07,-4.73466e-11,68841.1,21.2809], Tmin=(100,'K'), Tmax=(890.641,'K')), NASAPolynomial(coeffs=[4.29645,0.0226555,-1.15626e-05,2.16119e-09,-1.43083e-13,69565.2,15.3771], Tmin=(890.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=[C]OO[CH][C]([CH2])[O](25822)',
    structure = SMILES('[CH]=[C]OO[CH][C]([CH2])[O]'),
    E0 = (1015.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,250.644,1548.05],'cm^-1')),
        HinderedRotor(inertia=(0.13527,'amu*angstrom^2'), symmetry=1, barrier=(6.03075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135286,'amu*angstrom^2'), symmetry=1, barrier=(6.03076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0642,'amu*angstrom^2'), symmetry=1, barrier=(47.4367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135322,'amu*angstrom^2'), symmetry=1, barrier=(6.03075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06412,'amu*angstrom^2'), symmetry=1, barrier=(47.4368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.147143,0.103897,-0.000189456,1.72833e-07,-5.96746e-11,122314,36.5273], Tmin=(100,'K'), Tmax=(865.472,'K')), NASAPolynomial(coeffs=[9.8415,0.0307206,-1.58131e-05,3.01563e-09,-2.04378e-13,121597,-4.37448], Tmin=(865.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1015.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOOC) + radical(CJCO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=COO[CH][C]([CH2])[O](25719)',
    structure = SMILES('[CH]=COO[CH][C]([CH2])[O]'),
    E0 = (776.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,316.806,316.808],'cm^-1')),
        HinderedRotor(inertia=(0.620899,'amu*angstrom^2'), symmetry=1, barrier=(44.2219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24569,'amu*angstrom^2'), symmetry=1, barrier=(17.4986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.620902,'amu*angstrom^2'), symmetry=1, barrier=(44.2218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0790207,'amu*angstrom^2'), symmetry=1, barrier=(5.62797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0790197,'amu*angstrom^2'), symmetry=1, barrier=(5.62799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.312332,0.101572,-0.000160398,1.30957e-07,-4.20207e-11,93491.2,33.9944], Tmin=(100,'K'), Tmax=(816.046,'K')), NASAPolynomial(coeffs=[13.8173,0.0268863,-1.31394e-05,2.50571e-09,-1.72176e-13,91365.8,-30.1933], Tmin=(816.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(776.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOOC) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OOC[C]([CH2])[O](26857)',
    structure = SMILES('[CH]=[C]OOC[C]([CH2])[O]'),
    E0 = (829.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,246.318,2530.45],'cm^-1')),
        HinderedRotor(inertia=(0.20059,'amu*angstrom^2'), symmetry=1, barrier=(8.63632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200579,'amu*angstrom^2'), symmetry=1, barrier=(8.63632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86195,'amu*angstrom^2'), symmetry=1, barrier=(37.1118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200586,'amu*angstrom^2'), symmetry=1, barrier=(8.63631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861957,'amu*angstrom^2'), symmetry=1, barrier=(37.1117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.160744,0.102124,-0.00017559,1.56022e-07,-5.33378e-11,99892.9,36.2811], Tmin=(100,'K'), Tmax=(852.021,'K')), NASAPolynomial(coeffs=[10.456,0.0318801,-1.60088e-05,3.05494e-09,-2.08425e-13,98824.2,-8.89293], Tmin=(852.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(829.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OO[CH]C([CH2])[O](26858)',
    structure = SMILES('[CH]=[C]OO[CH]C([CH2])[O]'),
    E0 = (839.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,203.344,2292.5],'cm^-1')),
        HinderedRotor(inertia=(1.25768,'amu*angstrom^2'), symmetry=1, barrier=(36.9029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447605,'amu*angstrom^2'), symmetry=1, barrier=(13.1336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25768,'amu*angstrom^2'), symmetry=1, barrier=(36.9029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447605,'amu*angstrom^2'), symmetry=1, barrier=(13.1337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00989426,'amu*angstrom^2'), symmetry=1, barrier=(36.9029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0486047,0.0946556,-0.000150575,1.28228e-07,-4.30945e-11,101068,36.4017], Tmin=(100,'K'), Tmax=(818.83,'K')), NASAPolynomial(coeffs=[10.9792,0.0310614,-1.53962e-05,2.95909e-09,-2.04323e-13,99620.1,-12.0595], Tmin=(818.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(839.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCsJOOC) + radical(CJCO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OO[CH][C]([CH2])O(26859)',
    structure = SMILES('[CH]=[C]OO[CH][C]([CH2])O'),
    E0 = (785.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.307696,0.107859,-0.000195216,1.77612e-07,-6.11319e-11,94613.4,37.2957], Tmin=(100,'K'), Tmax=(869.721,'K')), NASAPolynomial(coeffs=[9.9258,0.0324804,-1.63797e-05,3.09853e-09,-2.0898e-13,93904.2,-4.49049], Tmin=(869.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(CCsJOOC) + radical(CJCO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OO[CH][C](C)[O](26860)',
    structure = SMILES('[CH]=[C]OO[CH][C](C)[O]'),
    E0 = (804.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,262.392,1871.14],'cm^-1')),
        HinderedRotor(inertia=(0.15578,'amu*angstrom^2'), symmetry=1, barrier=(7.61108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925051,'amu*angstrom^2'), symmetry=1, barrier=(45.1953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155781,'amu*angstrom^2'), symmetry=1, barrier=(7.61107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925043,'amu*angstrom^2'), symmetry=1, barrier=(45.1952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15578,'amu*angstrom^2'), symmetry=1, barrier=(7.61108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0253901,0.0984749,-0.000168942,1.52822e-07,-5.32511e-11,96861.1,35.2741], Tmin=(100,'K'), Tmax=(848.321,'K')), NASAPolynomial(coeffs=[8.89943,0.0347232,-1.74772e-05,3.34653e-09,-2.29104e-13,96143.8,-1.43181], Tmin=(848.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(804.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOOC) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (768.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (768.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1113.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1178.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1255.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1217.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (777.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (784.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (768.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1227.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (881.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (862.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (990.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (936.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (936.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    products = ['ketene(1375)', 'C=C([O])[CH][O](2850)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    products = ['[CH2]C([O])=COO[C]=C(25530)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]O[O](1354)', '[CH][C]([CH2])[O](13252)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H2CC(T)(1341)', '[CH2][C]([O])[CH]O[O](22833)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2][C]([O])[C]OO[C]=C(26854)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH][C]([O])[CH]OO[C]=C(26855)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    products = ['[CH2][C]([O])C1OOC1=C(26856)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]=[C]OOC=C([CH2])[O](25806)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=O(1376)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(36.0503,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination
Ea raised from 30.7 to 36.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]=[C]OO[CH][C]([CH2])[O](25822)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=COO[CH][C]([CH2])[O](25719)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OOC[C]([CH2])[O](26857)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]OO[CH]C([CH2])[O](26858)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]OO[CH][C]([CH2])O(26859)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]OO[CH][C](C)[O](26860)'],
    products = ['[CH2][C]([O])[CH]OO[C]=C(25548)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #6192',
    isomers = [
        '[CH2][C]([O])[CH]OO[C]=C(25548)',
    ],
    reactants = [
        ('ketene(1375)', 'C=C([O])[CH][O](2850)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #6192',
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


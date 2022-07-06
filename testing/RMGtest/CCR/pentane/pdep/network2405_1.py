species(
    label = '[CH]=COOC#C(2882)',
    structure = SMILES('[CH]=COOC#C'),
    E0 = (541.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.47296,'amu*angstrom^2'), symmetry=1, barrier=(33.8664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47269,'amu*angstrom^2'), symmetry=1, barrier=(33.8602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42267,0.0542774,-6.11083e-05,3.35115e-08,-7.15556e-12,65253.4,21.2931], Tmin=(100,'K'), Tmax=(1147.63,'K')), NASAPolynomial(coeffs=[13.7133,0.0114389,-5.11653e-06,9.85438e-10,-7.00658e-14,62432.3,-39.6948], Tmin=(1147.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'ketene(T)(1373)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (150.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32242,0.00993891,8.31e-06,-1.94096e-08,8.42706e-12,18150.6,9.4742], Tmin=(100,'K'), Tmax=(967.893,'K')), NASAPolynomial(coeffs=[7.72276,0.00396062,-1.3429e-06,2.69348e-10,-2.14899e-14,16727,-14.5651], Tmin=(967.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ketene(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[C]=COOC#C(10412)',
    structure = SMILES('[C]=COOC#C'),
    E0 = (852.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.4457,'amu*angstrom^2'), symmetry=1, barrier=(33.2396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44601,'amu*angstrom^2'), symmetry=1, barrier=(33.2465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44459,'amu*angstrom^2'), symmetry=1, barrier=(33.2139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60506,0.0553143,-7.48519e-05,5.02382e-08,-1.32664e-11,102648,20.621], Tmin=(100,'K'), Tmax=(927.699,'K')), NASAPolynomial(coeffs=[11.3583,0.0132614,-6.85694e-06,1.37582e-09,-9.8929e-14,100838,-25.7007], Tmin=(927.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(852.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'oxirene(6860)',
    structure = SMILES('C1=CO1'),
    E0 = (259.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31417,0.0119135,-5.00134e-07,-9.18335e-09,4.80348e-12,31257.7,8.09654], Tmin=(100,'K'), Tmax=(945.068,'K')), NASAPolynomial(coeffs=[6.90182,0.00459354,-1.36476e-06,2.32146e-10,-1.65654e-14,30228.3,-10.8673], Tmin=(945.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""oxirene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[C]1=CC=COO1(10413)',
    structure = SMILES('[C]1=CC=COO1'),
    E0 = (327.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6157,0.0162003,4.30322e-05,-7.02168e-08,2.8502e-11,39459.4,17.0487], Tmin=(100,'K'), Tmax=(959.896,'K')), NASAPolynomial(coeffs=[12.6799,0.0100005,-3.12765e-06,6.29672e-10,-5.16734e-14,35880.8,-39.6698], Tmin=(959.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1C=COO1(10414)',
    structure = SMILES('[CH]=C1C=COO1'),
    E0 = (341.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58495,0.0115418,6.8291e-05,-1.06054e-07,4.39368e-11,41201.1,17.251], Tmin=(100,'K'), Tmax=(931.334,'K')), NASAPolynomial(coeffs=[16.2692,0.00294775,1.31477e-06,-2.60048e-10,9.30814e-15,36476,-59.4775], Tmin=(931.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Cds_P)"""),
)

species(
    label = 'C#COOC#C(9631)',
    structure = SMILES('C#COOC#C'),
    E0 = (526.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200],'cm^-1')),
        HinderedRotor(inertia=(1.71762,'amu*angstrom^2'), symmetry=1, barrier=(39.4916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71352,'amu*angstrom^2'), symmetry=1, barrier=(39.3971,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63028,0.0552718,-7.33775e-05,4.90096e-08,-1.29944e-11,63367.3,17.7627], Tmin=(100,'K'), Tmax=(920.296,'K')), NASAPolynomial(coeffs=[10.8682,0.0151195,-7.9322e-06,1.60022e-09,-1.15391e-13,61667.1,-26.0374], Tmin=(920.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-OsCt) + group(Ct-CtOs) + group(Ct-CtOs) + group(Ct-CtH) + group(Ct-CtH)"""),
)

species(
    label = 'C#CO[O](2880)',
    structure = SMILES('C#CO[O]'),
    E0 = (341.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,207.446,207.513],'cm^-1')),
        HinderedRotor(inertia=(1.42756,'amu*angstrom^2'), symmetry=1, barrier=(43.6063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84093,0.0270252,-4.16758e-05,3.32096e-08,-1.02948e-11,41097.9,11.7896], Tmin=(100,'K'), Tmax=(872.981,'K')), NASAPolynomial(coeffs=[6.69652,0.00695047,-3.0441e-06,5.47527e-10,-3.61761e-14,40516.5,-5.76197], Tmin=(872.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-OsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08745,0.00578572,8.56332e-06,-1.72824e-08,7.83594e-12,26273.8,4.4608], Tmin=(100,'K'), Tmax=(897.945,'K')), NASAPolynomial(coeffs=[5.89068,0.00209055,4.88462e-08,-5.66903e-11,4.15085e-15,25415.9,-10.7352], Tmin=(897.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.303,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'HC2(2881)',
    structure = SMILES('[C]#C'),
    E0 = (556.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0139907,-3.08144e-05,3.10836e-08,-1.10946e-11,66983,5.75946], Tmin=(100,'K'), Tmax=(918.723,'K')), NASAPolynomial(coeffs=[3.14385,0.00498487,-2.32624e-06,4.08757e-10,-2.55787e-14,67315.6,7.08554], Tmin=(918.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.809,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""HC2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=CO[O](1355)',
    structure = SMILES('[CH]=CO[O]'),
    E0 = (345.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,3919.4],'cm^-1')),
        HinderedRotor(inertia=(0.735327,'amu*angstrom^2'), symmetry=1, barrier=(16.9066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89265,0.0254964,-3.40139e-05,2.50833e-08,-7.36731e-12,41533.8,14.7902], Tmin=(100,'K'), Tmax=(835.177,'K')), NASAPolynomial(coeffs=[6.45723,0.00842449,-3.35281e-06,6.09039e-10,-4.13537e-14,40938.3,-1.76488], Tmin=(835.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OOC#C(9626)',
    structure = SMILES('[CH]=[C]OOC#C'),
    E0 = (781.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,2175,525,3120,650,792.5,1650,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.39565,'amu*angstrom^2'), symmetry=1, barrier=(32.0887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39392,'amu*angstrom^2'), symmetry=1, barrier=(32.0491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78262,0.053915,-7.85851e-05,5.62666e-08,-1.41104e-11,94068,23.1506], Tmin=(100,'K'), Tmax=(631.548,'K')), NASAPolynomial(coeffs=[8.84699,0.0168337,-8.71075e-06,1.71653e-09,-1.20853e-13,93023,-8.89392], Tmin=(631.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(781.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]#COOC=[CH](10415)',
    structure = SMILES('[C]#COOC=[CH]'),
    E0 = (878.902,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.48697,'amu*angstrom^2'), symmetry=1, barrier=(34.1883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48913,'amu*angstrom^2'), symmetry=1, barrier=(34.2381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48797,'amu*angstrom^2'), symmetry=1, barrier=(34.2115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46499,0.0575612,-8.26829e-05,5.90679e-08,-1.64085e-11,105797,21.8859], Tmin=(100,'K'), Tmax=(888.758,'K')), NASAPolynomial(coeffs=[11.6642,0.0116607,-5.21916e-06,9.65269e-10,-6.57398e-14,103984,-26.1174], Tmin=(888.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(878.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_P) + radical(Acetyl)"""),
)

species(
    label = 'C#COO[C]=C(2664)',
    structure = SMILES('C#COO[C]=C'),
    E0 = (534.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.34249,'amu*angstrom^2'), symmetry=1, barrier=(30.8664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35112,'amu*angstrom^2'), symmetry=1, barrier=(31.0648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34311,'amu*angstrom^2'), symmetry=1, barrier=(30.8808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88133,0.0503015,-6.29738e-05,4.33221e-08,-1.22468e-11,64347.3,22.3829], Tmin=(100,'K'), Tmax=(853.703,'K')), NASAPolynomial(coeffs=[8.31092,0.0201755,-1.00404e-05,1.98546e-09,-1.41558e-13,63249.5,-7.61908], Tmin=(853.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = '[C]#COOC=C(2883)',
    structure = SMILES('[C]#COOC=C'),
    E0 = (631.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2175,525,180],'cm^-1')),
        HinderedRotor(inertia=(1.48049,'amu*angstrom^2'), symmetry=1, barrier=(34.0394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47997,'amu*angstrom^2'), symmetry=1, barrier=(34.0275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47954,'amu*angstrom^2'), symmetry=1, barrier=(34.0175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57405,0.0536279,-6.46079e-05,3.98496e-08,-9.69897e-12,76075.9,21.091], Tmin=(100,'K'), Tmax=(1005.47,'K')), NASAPolynomial(coeffs=[11.4043,0.0145206,-6.26589e-06,1.1665e-09,-8.07773e-14,74099,-26.3881], Tmin=(1005.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
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
    E0 = (577.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1064.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (618.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (564.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (557.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (753.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (592.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (541.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (932.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (901.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (993.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1090.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (647.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (694.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=COOC#C(2882)'],
    products = ['HCCO(2227)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.29914e+17,'s^-1'), n=-1.73308, Ea=(36.1101,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.16692063000437535, var=8.814681701569807, Tref=1000.0, N=65, correlation='Root',), comment="""BM rule fitted to 2 training reactions at node Root
    Total Standard Deviation in ln(k): 6.371362717619001
Exact match found for rate rule [Root]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[C]=COOC#C(10412)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=COOC#C(2882)'],
    products = ['HCCO(2227)', 'oxirene(6860)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.88184e+13,'s^-1'), n=-0.182188, Ea=(76.9763,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO;Y_rad_intra;OOR] for rate rule [R2OO_D;Cd_pri_rad_in;OOR]
Euclidian distance = 2.23606797749979
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=COOC#C(2882)'],
    products = ['[C]1=CC=COO1(10413)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.926e+10,'s^-1'), n=0.198, Ea=(22.8237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_DSR;triplebond_intra_H;radadd_intra_cdsingleH] for rate rule [R6_DSS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=COOC#C(2882)'],
    products = ['[CH]=C1C=COO1(10414)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.42e+11,'s^-1'), n=0.258, Ea=(15.8866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_DSR;triplebond_intra_H;radadd_intra_cdsingleH] for rate rule [R6_DSS_T;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#COOC#C(9631)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8556.54,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C#CO[O](2880)', 'C2H2(1342)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.04e+12,'cm^3/(mol*s)'), n=0, Ea=(33.0536,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Ct-H_Ct-H;O_rad/NonDe] for rate rule [Ct-H_Ct-H;OJ-O2s]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCCO(2227)', 'ketene(T)(1373)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(225.895,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
family: R_Recombination
Ea raised from 221.6 to 225.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CO[O](2880)', 'C2H2(T)(1343)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(212954,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['HC2(2881)', '[CH]=CO[O](1355)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]=[C]OOC#C(9626)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.88213e+06,'m^3/(mol*s)'), n=0.314663, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.379272271586, var=0.88677526262, Tref=1000.0, N=15, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O
    Total Standard Deviation in ln(k): 2.84077927867
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[C]#COOC=[CH](10415)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=COOC#C(2882)'],
    products = ['C#COO[C]=C(2664)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]#COOC=C(2883)'],
    products = ['[CH]=COOC#C(2882)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(26171.5,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_TSSSD;Ct_rad_out;Cd_H_out_singleH]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2405',
    isomers = [
        '[CH]=COOC#C(2882)',
    ],
    reactants = [
        ('HCCO(2227)', 'ketene(1375)'),
        ('HCCO(2227)', 'ketene(T)(1373)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2405',
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


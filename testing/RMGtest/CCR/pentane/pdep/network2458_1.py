species(
    label = '[CH2][CH]C(C)[CH][C]=O(4178)',
    structure = SMILES('[CH2][CH]C(C)[CH][C]=O'),
    E0 = (446.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,850.486,850.567],'cm^-1')),
        HinderedRotor(inertia=(0.122248,'amu*angstrom^2'), symmetry=1, barrier=(2.81072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122387,'amu*angstrom^2'), symmetry=1, barrier=(2.81391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00547187,'amu*angstrom^2'), symmetry=1, barrier=(2.80925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122372,'amu*angstrom^2'), symmetry=1, barrier=(2.81356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122311,'amu*angstrom^2'), symmetry=1, barrier=(2.81218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37781,0.0597816,-4.89364e-05,2.32743e-08,-4.82459e-12,53844.7,30.1838], Tmin=(100,'K'), Tmax=(1105.4,'K')), NASAPolynomial(coeffs=[8.29099,0.0347652,-1.49891e-05,2.80038e-09,-1.94078e-13,52316.4,-3.86095], Tmin=(1105.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    label = 'm1_allyl(186)',
    structure = SMILES('C=C[CH]C'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800698,'amu*angstrom^2'), symmetry=1, barrier=(1.84096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.062381,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C(C)C=C=O(4176)',
    structure = SMILES('[CH2][CH]C(C)C=C=O'),
    E0 = (241.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,344.134,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0106282,'amu*angstrom^2'), symmetry=1, barrier=(26.2918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919169,'amu*angstrom^2'), symmetry=1, barrier=(7.72529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919256,'amu*angstrom^2'), symmetry=1, barrier=(7.72529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312855,'amu*angstrom^2'), symmetry=1, barrier=(26.2917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4928,0.058788,-4.53704e-05,1.9742e-08,-3.79184e-12,29088.5,28.2399], Tmin=(100,'K'), Tmax=(1162.17,'K')), NASAPolynomial(coeffs=[8.00864,0.0363616,-1.64249e-05,3.13769e-09,-2.20012e-13,27574,-4.17453], Tmin=(1162.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'CH2(S)(11)',
    structure = SMILES('[CH2]'),
    E0 = (418.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1381.44,2681.05,3148.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08586,-0.00146571,5.70225e-06,-3.92867e-09,8.88752e-13,50365.3,-0.343718], Tmin=(100,'K'), Tmax=(1285.36,'K')), NASAPolynomial(coeffs=[2.65113,0.0038016,-1.38109e-06,2.30893e-10,-1.47438e-14,50667.9,6.68035], Tmin=(1285.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C[CH][C]=O(2616)',
    structure = SMILES('[CH2][CH]C[CH][C]=O'),
    E0 = (476.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,1652.74,1652.79],'cm^-1')),
        HinderedRotor(inertia=(0.1661,'amu*angstrom^2'), symmetry=1, barrier=(8.24258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333532,'amu*angstrom^2'), symmetry=1, barrier=(64.6535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00425184,'amu*angstrom^2'), symmetry=1, barrier=(8.24363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16632,'amu*angstrom^2'), symmetry=1, barrier=(8.24369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3428.95,'J/mol'), sigma=(5.98513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.59 K, Pc=36.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8346,0.0511755,-6.13261e-05,4.91033e-08,-1.65214e-11,57360.1,26.7553], Tmin=(100,'K'), Tmax=(839.036,'K')), NASAPolynomial(coeffs=[5.14519,0.0296716,-1.2654e-05,2.30336e-09,-1.55372e-13,57005.9,12.5647], Tmin=(839.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(C)[CH][CH][C]=O(10257)',
    structure = SMILES('[CH2]C(C)[CH]C=[C][O]'),
    E0 = (445.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,552.602,552.824,3568.08],'cm^-1')),
        HinderedRotor(inertia=(0.0916064,'amu*angstrom^2'), symmetry=1, barrier=(19.9569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.384347,'amu*angstrom^2'), symmetry=1, barrier=(83.404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0110255,'amu*angstrom^2'), symmetry=1, barrier=(2.41842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.62926,'amu*angstrom^2'), symmetry=1, barrier=(83.4438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04334,0.0558583,-2.37753e-05,-9.87339e-09,8.31093e-12,53694.1,29.9723], Tmin=(100,'K'), Tmax=(977.435,'K')), NASAPolynomial(coeffs=[13.3851,0.0258418,-9.15628e-06,1.603e-09,-1.10043e-13,50302.7,-34.2951], Tmin=(977.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][CH]C(C)[C]=O(10258)',
    structure = SMILES('[CH2][CH][CH]C(C)[C]=O'),
    E0 = (477.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,720.454,2237.6],'cm^-1')),
        HinderedRotor(inertia=(0.00330972,'amu*angstrom^2'), symmetry=1, barrier=(11.7595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00330964,'amu*angstrom^2'), symmetry=1, barrier=(11.7592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0102639,'amu*angstrom^2'), symmetry=1, barrier=(3.7803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.45309,'amu*angstrom^2'), symmetry=1, barrier=(56.4013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164415,'amu*angstrom^2'), symmetry=1, barrier=(3.78022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50167,0.0588729,-5.32842e-05,3.18483e-08,-8.71876e-12,57559.7,33.1439], Tmin=(100,'K'), Tmax=(841.978,'K')), NASAPolynomial(coeffs=[5.76945,0.0385976,-1.71629e-05,3.24751e-09,-2.26501e-13,56841,13.2884], Tmin=(841.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
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
    label = 'C[CH][CH][C]=O(9638)',
    structure = SMILES('C[CH]C=[C][O]'),
    E0 = (297.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,484.767,484.8],'cm^-1')),
        HinderedRotor(inertia=(0.000717709,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163781,'amu*angstrom^2'), symmetry=1, barrier=(27.3151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43094,0.0256217,1.17561e-05,-3.33743e-08,1.45333e-11,35785.9,19.0015], Tmin=(100,'K'), Tmax=(974.935,'K')), NASAPolynomial(coeffs=[10.1488,0.0154088,-5.53624e-06,1.01957e-09,-7.35997e-14,33261.5,-23.2652], Tmin=(974.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C][O](6861)',
    structure = SMILES('[CH][C]=O'),
    E0 = (425.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1491.4,1492.35],'cm^-1')),
        HinderedRotor(inertia=(0.0736325,'amu*angstrom^2'), symmetry=1, barrier=(5.08901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62981,0.00805288,-4.35162e-06,1.03741e-09,-8.43183e-14,51134.2,10.4631], Tmin=(100,'K'), Tmax=(1941.7,'K')), NASAPolynomial(coeffs=[6.62042,0.00292938,-1.19496e-06,2.28726e-10,-1.56225e-14,49777.3,-6.4529], Tmin=(1941.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH][CH]C(3856)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (450.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2810.32],'cm^-1')),
        HinderedRotor(inertia=(0.00169393,'amu*angstrom^2'), symmetry=1, barrier=(9.48541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010072,'amu*angstrom^2'), symmetry=1, barrier=(56.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100766,'amu*angstrom^2'), symmetry=1, barrier=(56.4312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79546,0.0335488,-4.51222e-05,5.0661e-08,-2.13657e-11,54172.4,20.3345], Tmin=(100,'K'), Tmax=(854.972,'K')), NASAPolynomial(coeffs=[-1.26875,0.0345785,-1.53757e-05,2.86235e-09,-1.94746e-13,55524.7,43.1492], Tmin=(854.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][C]C(C)[CH][C]=O(10259)',
    structure = SMILES('[CH2][C]C(C)[CH][C]=O'),
    E0 = (700.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961581,0.0661273,-6.44529e-05,3.45686e-08,-7.5847e-12,84371.9,29.0245], Tmin=(100,'K'), Tmax=(1092.92,'K')), NASAPolynomial(coeffs=[11.5872,0.0272383,-1.10788e-05,2.01105e-09,-1.37317e-13,82049.4,-23.1821], Tmin=(1092.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)[C][C]=O(10260)',
    structure = SMILES('[CH2][CH]C(C)[C][C]=O'),
    E0 = (727.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00387,0.0710864,-9.28035e-05,7.51447e-08,-2.55208e-11,87613.7,31.669], Tmin=(100,'K'), Tmax=(770.223,'K')), NASAPolynomial(coeffs=[7.18821,0.0354372,-1.64982e-05,3.14478e-09,-2.18401e-13,86765.9,4.1278], Tmin=(770.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][CH]C(C)[CH][C]=O(10261)',
    structure = SMILES('[CH][CH]C(C)[CH][C]=O'),
    E0 = (689.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,207.429,804.968,1073.29,1382.86,1626.99],'cm^-1')),
        HinderedRotor(inertia=(0.147192,'amu*angstrom^2'), symmetry=1, barrier=(3.55446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147192,'amu*angstrom^2'), symmetry=1, barrier=(3.55446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147192,'amu*angstrom^2'), symmetry=1, barrier=(3.55446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147192,'amu*angstrom^2'), symmetry=1, barrier=(3.55446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147192,'amu*angstrom^2'), symmetry=1, barrier=(3.55446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24713,0.0625527,-5.94395e-05,3.24433e-08,-7.48231e-12,83072,29.8729], Tmin=(100,'K'), Tmax=(1023.57,'K')), NASAPolynomial(coeffs=[9.26812,0.0312078,-1.35053e-05,2.52585e-09,-1.75232e-13,81429.9,-9.01061], Tmin=(1023.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1C(C)C1[C]=O(10262)',
    structure = SMILES('[CH2]C1C(C)C1[C]=O'),
    E0 = (227.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45485,0.0421531,1.67831e-05,-5.3841e-08,2.45052e-11,27509.4,25.668], Tmin=(100,'K'), Tmax=(948.151,'K')), NASAPolynomial(coeffs=[14.1171,0.0230529,-7.29272e-06,1.26189e-09,-8.94119e-14,23565.7,-42.8806], Tmin=(948.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]CC(C)=C[C]=O(10263)',
    structure = SMILES('[CH2]CC(C)=C[C]=O'),
    E0 = (186.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909237,0.0744654,-8.64557e-05,6.11295e-08,-1.86581e-11,22572.7,26.3375], Tmin=(100,'K'), Tmax=(778.83,'K')), NASAPolynomial(coeffs=[7.67038,0.0397378,-1.9566e-05,3.86814e-09,-2.76007e-13,21519.6,-4.59046], Tmin=(778.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C=C(C)C[C]=O(4175)',
    structure = SMILES('C=C[C](C)C[C]=O'),
    E0 = (137.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,2367.61],'cm^-1')),
        HinderedRotor(inertia=(0.00423044,'amu*angstrom^2'), symmetry=1, barrier=(16.8281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164387,'amu*angstrom^2'), symmetry=1, barrier=(16.8269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.468175,'amu*angstrom^2'), symmetry=1, barrier=(10.7643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.50028,'amu*angstrom^2'), symmetry=1, barrier=(80.4783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31246,0.0562455,-3.87501e-05,1.43709e-08,-2.24778e-12,16586.7,25.0716], Tmin=(100,'K'), Tmax=(1449.24,'K')), NASAPolynomial(coeffs=[10.5775,0.0306736,-1.22827e-05,2.19568e-09,-1.47508e-13,13901.2,-23.0644], Tmin=(1449.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C](C)[CH][C]=O(9983)',
    structure = SMILES('C=C[C](C)[CH][C]=O'),
    E0 = (304.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,751.629,753.362],'cm^-1')),
        HinderedRotor(inertia=(3.44203,'amu*angstrom^2'), symmetry=1, barrier=(79.139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0023116,'amu*angstrom^2'), symmetry=1, barrier=(13.6128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016343,'amu*angstrom^2'), symmetry=1, barrier=(5.0883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.44192,'amu*angstrom^2'), symmetry=1, barrier=(79.1364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34258,0.0514931,-2.72118e-05,2.73287e-10,3.40644e-12,36738.7,24.7245], Tmin=(100,'K'), Tmax=(1018,'K')), NASAPolynomial(coeffs=[11.1854,0.0262021,-9.66737e-06,1.69886e-09,-1.15331e-13,34041.2,-26.3433], Tmin=(1018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C[C]=O(9636)',
    structure = SMILES('C[CH]C=C=O'),
    E0 = (39.4293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,1069.02],'cm^-1')),
        HinderedRotor(inertia=(0.163668,'amu*angstrom^2'), symmetry=1, barrier=(3.76304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30969,'amu*angstrom^2'), symmetry=1, barrier=(30.1123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47551,0.0278456,-3.33541e-06,-1.12635e-08,5.08824e-12,4801.91,15.1824], Tmin=(100,'K'), Tmax=(1090.48,'K')), NASAPolynomial(coeffs=[7.95588,0.0195027,-8.03539e-06,1.49913e-09,-1.04907e-13,2907.47,-14.9378], Tmin=(1090.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.4293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_S)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([532.913,1391.12,1391.12,2779.21,3448.45,3448.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.948,0.0008276,8.34932e-06,-9.82634e-09,3.80104e-12,16425.4,0.336655], Tmin=(100,'K'), Tmax=(660.467,'K')), NASAPolynomial(coeffs=[3.2217,0.00522646,-1.64125e-06,2.58225e-10,-1.62579e-14,16521.3,3.53938], Tmin=(660.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C=C[CH][C]=O(9826)',
    structure = SMILES('[CH2]C=C[CH][C]=O'),
    E0 = (344.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,579.211],'cm^-1')),
        HinderedRotor(inertia=(0.109865,'amu*angstrom^2'), symmetry=1, barrier=(26.367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109836,'amu*angstrom^2'), symmetry=1, barrier=(26.387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110789,'amu*angstrom^2'), symmetry=1, barrier=(26.3527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96036,0.0299846,2.4079e-05,-5.41179e-08,2.23924e-11,41511.2,22.2791], Tmin=(100,'K'), Tmax=(1008.12,'K')), NASAPolynomial(coeffs=[14.799,0.0148331,-6.62881e-06,1.40459e-09,-1.09437e-13,37104,-48.7837], Tmin=(1008.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH][CH][C]=O(9827)',
    structure = SMILES('[CH2][CH][CH]C=[C][O]'),
    E0 = (672.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,334.418,334.876,3283.27],'cm^-1')),
        HinderedRotor(inertia=(0.00368016,'amu*angstrom^2'), symmetry=1, barrier=(3.7932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475898,'amu*angstrom^2'), symmetry=1, barrier=(49.0244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0948626,'amu*angstrom^2'), symmetry=1, barrier=(97.7017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11421,0.0402686,-2.82926e-05,1.07211e-08,-1.72234e-12,81002,27.5963], Tmin=(100,'K'), Tmax=(1406.92,'K')), NASAPolynomial(coeffs=[8.42123,0.0223372,-9.17472e-06,1.66212e-09,-1.12619e-13,79227.3,-4.98465], Tmin=(1406.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][C](C)[CH][C]=O(10264)',
    structure = SMILES('[CH2][CH][C](C)C=[C][O]'),
    E0 = (632.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,365.404,365.497,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0383423,'amu*angstrom^2'), symmetry=1, barrier=(47.7038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0632,'amu*angstrom^2'), symmetry=1, barrier=(100.789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252662,'amu*angstrom^2'), symmetry=1, barrier=(28.6874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010965,'amu*angstrom^2'), symmetry=1, barrier=(28.6871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23285,0.0542361,-3.77105e-05,1.35982e-08,-2.00723e-12,76132.6,31.9518], Tmin=(100,'K'), Tmax=(1570.11,'K')), NASAPolynomial(coeffs=[12.7485,0.0248987,-9.68298e-06,1.69767e-09,-1.12373e-13,72516.4,-28.8], Tmin=(1570.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Cs_S) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][C]=O(9158)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][C]=O'),
    E0 = (652.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,944.543,944.577],'cm^-1')),
        HinderedRotor(inertia=(0.00506013,'amu*angstrom^2'), symmetry=1, barrier=(3.2039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139316,'amu*angstrom^2'), symmetry=1, barrier=(3.20315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00506097,'amu*angstrom^2'), symmetry=1, barrier=(3.20391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139369,'amu*angstrom^2'), symmetry=1, barrier=(3.20436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139367,'amu*angstrom^2'), symmetry=1, barrier=(3.20431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20527,0.0636919,-6.9144e-05,4.57284e-08,-1.26958e-11,78516.4,32.586], Tmin=(100,'K'), Tmax=(865.982,'K')), NASAPolynomial(coeffs=[8.25811,0.0311141,-1.27137e-05,2.28552e-09,-1.54076e-13,77294.9,-0.42492], Tmin=(865.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C[C](C)[CH][C]=O(10265)',
    structure = SMILES('[CH2]C[C](C)C=[C][O]'),
    E0 = (437.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11239,0.0547871,-2.47862e-05,-5.36683e-09,5.68842e-12,52741.1,29.1247], Tmin=(100,'K'), Tmax=(1024.08,'K')), NASAPolynomial(coeffs=[12.5356,0.0275714,-1.04127e-05,1.86973e-09,-1.29016e-13,49488.9,-30.7127], Tmin=(1024.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][C](C)C[C]=O(4177)',
    structure = SMILES('[CH2][CH][C](C)C[C]=O'),
    E0 = (464.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,2648.03,2648.04],'cm^-1')),
        HinderedRotor(inertia=(0.145705,'amu*angstrom^2'), symmetry=1, barrier=(3.35005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145676,'amu*angstrom^2'), symmetry=1, barrier=(3.34938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145659,'amu*angstrom^2'), symmetry=1, barrier=(3.34898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145719,'amu*angstrom^2'), symmetry=1, barrier=(3.35037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145656,'amu*angstrom^2'), symmetry=1, barrier=(3.34891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08089,0.0762624,-0.000123425,1.2233e-07,-4.65259e-11,55997.6,32.6163], Tmin=(100,'K'), Tmax=(848.661,'K')), NASAPolynomial(coeffs=[1.27266,0.0473271,-2.27369e-05,4.31428e-09,-2.95271e-13,56974.5,37.6699], Tmin=(848.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC([CH2])[CH][C]=O(5300)',
    structure = SMILES('[CH2]CC([CH2])[CH][C]=O'),
    E0 = (457.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,338.159,338.16],'cm^-1')),
        HinderedRotor(inertia=(0.00147421,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00147422,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00147419,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0014742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00147422,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922897,0.0660645,-6.20976e-05,3.35322e-08,-7.47196e-12,55132,31.0378], Tmin=(100,'K'), Tmax=(1076.22,'K')), NASAPolynomial(coeffs=[10.7763,0.0294415,-1.10529e-05,1.91188e-09,-1.26599e-13,53011.1,-17.2228], Tmin=(1076.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[C]=O(4179)',
    structure = SMILES('[CH2][CH]C([CH2])C[C]=O'),
    E0 = (484.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,180,1543.7],'cm^-1')),
        HinderedRotor(inertia=(0.135513,'amu*angstrom^2'), symmetry=1, barrier=(3.1157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135936,'amu*angstrom^2'), symmetry=1, barrier=(3.12545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136415,'amu*angstrom^2'), symmetry=1, barrier=(3.13646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00185875,'amu*angstrom^2'), symmetry=1, barrier=(3.13707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136822,'amu*angstrom^2'), symmetry=1, barrier=(3.14581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01958,0.0703602,-8.79602e-05,7.04702e-08,-2.35795e-11,58371.4,33.4881], Tmin=(100,'K'), Tmax=(820.932,'K')), NASAPolynomial(coeffs=[6.49036,0.037436,-1.6349e-05,3.01548e-09,-2.05119e-13,57684.4,9.46055], Tmin=(820.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][C](C)[CH][C]=O(10149)',
    structure = SMILES('C[CH][C](C)C=[C][O]'),
    E0 = (426.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,349.54,349.653,349.662],'cm^-1')),
        HinderedRotor(inertia=(4.10431e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00137927,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00137888,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46318,'amu*angstrom^2'), symmetry=1, barrier=(40.1616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21737,0.0533312,-2.70727e-05,1.70212e-09,2.02869e-12,51448.9,29.5223], Tmin=(100,'K'), Tmax=(1134.45,'K')), NASAPolynomial(coeffs=[11.4077,0.0295571,-1.17118e-05,2.12105e-09,-1.45249e-13,48354.6,-24.3734], Tmin=(1134.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Cs_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][C](C)[CH]C=O(10266)',
    structure = SMILES('[CH2][CH][C](C)C=C[O]'),
    E0 = (392.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12027,0.0515094,-8.34458e-06,-2.6998e-08,1.43422e-11,47307.2,29.2164], Tmin=(100,'K'), Tmax=(982.821,'K')), NASAPolynomial(coeffs=[14.4169,0.0246471,-8.94204e-06,1.62215e-09,-1.14859e-13,43377.3,-41.398], Tmin=(982.821,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][C]=O(4240)',
    structure = SMILES('[CH2]C([CH]C)[CH][C]=O'),
    E0 = (446.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,957.83,958.279],'cm^-1')),
        HinderedRotor(inertia=(0.0744504,'amu*angstrom^2'), symmetry=1, barrier=(1.71176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749299,'amu*angstrom^2'), symmetry=1, barrier=(1.72279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749341,'amu*angstrom^2'), symmetry=1, barrier=(1.72288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0744698,'amu*angstrom^2'), symmetry=1, barrier=(1.71221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748947,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2049,0.0625364,-5.72642e-05,3.16551e-08,-7.47846e-12,53832.1,30.8], Tmin=(100,'K'), Tmax=(999.296,'K')), NASAPolynomial(coeffs=[8.4989,0.0333397,-1.3438e-05,2.41688e-09,-1.63698e-13,52374.4,-4.38418], Tmin=(999.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH]C=O(5533)',
    structure = SMILES('[CH2][CH]C([CH2])C=C[O]'),
    E0 = (465.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,726.327,3368.12],'cm^-1')),
        HinderedRotor(inertia=(0.0081838,'amu*angstrom^2'), symmetry=1, barrier=(3.08188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8863,'amu*angstrom^2'), symmetry=1, barrier=(20.3778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0547218,'amu*angstrom^2'), symmetry=1, barrier=(20.377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219919,'amu*angstrom^2'), symmetry=1, barrier=(82.3948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.872014,0.0584622,-2.5888e-05,-1.2245e-08,1.02813e-11,56107,32.9424], Tmin=(100,'K'), Tmax=(960.31,'K')), NASAPolynomial(coeffs=[15.3348,0.0224657,-7.53297e-06,1.30371e-09,-9.02495e-14,52211.3,-42.0672], Tmin=(960.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1C([O])=CC1C(10267)',
    structure = SMILES('[CH2]C1C(=O)[CH]C1C'),
    E0 = (194.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4889,0.0380655,3.49868e-05,-7.25304e-08,3.02453e-11,23480.4,23.8923], Tmin=(100,'K'), Tmax=(974.693,'K')), NASAPolynomial(coeffs=[15.0022,0.0245944,-8.89536e-06,1.67809e-09,-1.239e-13,18851.8,-51.1858], Tmin=(974.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C=C(C)C=C[O](10268)',
    structure = SMILES('[CH2]C=C(C)C=C[O]'),
    E0 = (71.0307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839611,0.0524779,8.11684e-06,-5.76576e-08,2.88436e-11,8672.47,24.5727], Tmin=(100,'K'), Tmax=(933.682,'K')), NASAPolynomial(coeffs=[19.2511,0.0168367,-4.08365e-06,6.49357e-10,-4.81251e-14,3349.84,-73.0805], Tmin=(933.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.0307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC(C)[C]=[C][O](9997)',
    structure = SMILES('C=CC(C)[C][C]=O'),
    E0 = (453.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,288.98,289.036,289.071],'cm^-1')),
        HinderedRotor(inertia=(0.00201806,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203706,'amu*angstrom^2'), symmetry=1, barrier=(12.0635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339522,'amu*angstrom^2'), symmetry=1, barrier=(20.114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203684,'amu*angstrom^2'), symmetry=1, barrier=(12.0639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14321,0.062784,-5.70781e-05,2.8431e-08,-5.86876e-12,54625,27.3343], Tmin=(100,'K'), Tmax=(1145.79,'K')), NASAPolynomial(coeffs=[10.9903,0.0284074,-1.20745e-05,2.24621e-09,-1.55523e-13,52368.4,-21.5126], Tmin=(1145.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)[C]=C[O](10269)',
    structure = SMILES('[CH2][CH]C(C)[C]=C[O]'),
    E0 = (498.246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180,652.674,4000],'cm^-1')),
        HinderedRotor(inertia=(0.438051,'amu*angstrom^2'), symmetry=1, barrier=(21.8817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.072362,'amu*angstrom^2'), symmetry=1, barrier=(21.8813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951707,'amu*angstrom^2'), symmetry=1, barrier=(21.8816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00757987,'amu*angstrom^2'), symmetry=1, barrier=(2.29103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76993,0.0619499,-4.0028e-05,6.69379e-09,2.10982e-12,60049.3,32.0025], Tmin=(100,'K'), Tmax=(1053.14,'K')), NASAPolynomial(coeffs=[14.7908,0.0246156,-9.52621e-06,1.73848e-09,-1.21073e-13,56213.3,-40.5574], Tmin=(1053.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(C)[C]=[C][O](10270)',
    structure = SMILES('[CH2]CC(C)[C][C]=O'),
    E0 = (533.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969303,0.0703013,-7.35263e-05,4.51688e-08,-1.17321e-11,64218.8,29.2469], Tmin=(100,'K'), Tmax=(916.093,'K')), NASAPolynomial(coeffs=[9.08886,0.0348485,-1.54766e-05,2.92463e-09,-2.03811e-13,62731.1,-9.21378], Tmin=(916.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)[C]=[C]O(10271)',
    structure = SMILES('[CH2][CH]C(C)[C]=[C]O'),
    E0 = (596.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1670,1700,300,440,339.882,339.892],'cm^-1')),
        HinderedRotor(inertia=(0.00145915,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00145938,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00145938,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147532,'amu*angstrom^2'), symmetry=1, barrier=(12.0939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147526,'amu*angstrom^2'), symmetry=1, barrier=(12.0938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691768,0.0692395,-6.53535e-05,3.29012e-08,-6.67506e-12,71867.8,34.0005], Tmin=(100,'K'), Tmax=(1187.1,'K')), NASAPolynomial(coeffs=[13.7077,0.0253814,-9.93519e-06,1.77867e-09,-1.20747e-13,68777.5,-31.0265], Tmin=(1187.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]C(C)[C][C]=O(10146)',
    structure = SMILES('C[CH]C(C)[C][C]=O'),
    E0 = (522.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10512,0.0685392,-7.50787e-05,5.19402e-08,-1.56318e-11,62925.2,29.5305], Tmin=(100,'K'), Tmax=(788.628,'K')), NASAPolynomial(coeffs=[7.05989,0.0383359,-1.76304e-05,3.37605e-09,-2.36536e-13,61986,2.21625], Tmin=(788.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][C](C)C=[C]O(10272)',
    structure = SMILES('[CH2]C=C(C)[CH][C]O'),
    E0 = (477.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202123,0.0806235,-8.14128e-05,4.2186e-08,-8.7129e-12,57626.3,25.1896], Tmin=(100,'K'), Tmax=(1171.04,'K')), NASAPolynomial(coeffs=[16.3202,0.025568,-1.08918e-05,2.03899e-09,-1.42134e-13,53851.3,-55.1159], Tmin=(1171.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])C=[C]O(10273)',
    structure = SMILES('[CH2][CH]C([CH2])C=[C]O'),
    E0 = (563.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,254.028,3075.6],'cm^-1')),
        HinderedRotor(inertia=(0.00261235,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00261233,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0026124,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474045,'amu*angstrom^2'), symmetry=1, barrier=(21.7064,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77016,'amu*angstrom^2'), symmetry=1, barrier=(81.0466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.662544,0.0672867,-5.65304e-05,2.08115e-08,-1.41523e-12,67931.3,35.4125], Tmin=(100,'K'), Tmax=(958.985,'K')), NASAPolynomial(coeffs=[14.5074,0.0228027,-7.69756e-06,1.28673e-09,-8.52186e-14,64666,-33.9812], Tmin=(958.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(C=CJO)"""),
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
    E0 = (446.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (446.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (895.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (644.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (677.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (910.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (930.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (912.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (939.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (901.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (449.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (469.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (469.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (524.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (546.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (597.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (514.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (809.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (843.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (863.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (605.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (623.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (575.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (602.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (604.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (600.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (563.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (596.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (455.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (510.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (680.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (637.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (702.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (675.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (746.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (691.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (561.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (637.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['HCCO(2227)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2][CH]C(C)C=C=O(4176)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]C[CH][C]=O(2616)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C)[CH][CH][C]=O(10257)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.37952e+08,'s^-1'), n=1.37167, Ea=(199.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][CH]C(C)[C]=O(10258)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.37952e+08,'s^-1'), n=1.37167, Ea=(199.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH2](502)', 'C[CH][CH][C]=O(9638)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C]C(C)[CH][C]=O(10259)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][CH]C(C)[C][C]=O(10260)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH][CH]C(C)[CH][C]=O(10261)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]C1C(C)C1[C]=O(10262)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]CC(C)=C[C]=O(10263)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]C=C(C)C[C]=O(4175)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'C=C[C](C)[CH][C]=O(9983)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C][O](6861)', 'm1_allyl(186)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH][CH2](502)', 'CC=C[C]=O(9636)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH3(17)', '[CH2]C=C[CH][C]=O(9826)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.139979,'m^3/(mol*s)'), n=2.09962, Ea=(33.817,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', '[CH2][CH][CH][CH][C]=O(9827)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2][CH][C](C)[CH][C]=O(10264)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[CH][C]=O(9158)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]C[C](C)[CH][C]=O(10265)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][C](C)C[C]=O(4177)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CC([CH2])[CH][C]=O(5300)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([CH2])C[C]=O(4179)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['C[CH][C](C)[CH][C]=O(10149)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2][CH][C](C)[CH]C=O(10266)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]C([CH]C)[CH][C]=O(4240)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(301993,'s^-1'), n=1.775, Ea=(116.315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R4HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]C([CH2])[CH]C=O(5533)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_2;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]C1C([O])=CC1C(10267)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['[CH2]C=C(C)C=C[O](10268)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', 'C=CC(C)[C]=[C][O](9997)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HCCO(2227)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C(C)[C]=C[O](10269)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]CC(C)[C]=[C][O](10270)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]C(C)[C]=[C]O(10271)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C[CH]C(C)[C][C]=O(10146)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(115.297,'s^-1'), n=2.99825, Ea=(169.04,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH][C](C)C=[C]O(10272)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C=[C]O(10273)'],
    products = ['[CH2][CH]C(C)[CH][C]=O(4178)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2458',
    isomers = [
        '[CH2][CH]C(C)[CH][C]=O(4178)',
    ],
    reactants = [
        ('HCCO(2227)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2458',
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


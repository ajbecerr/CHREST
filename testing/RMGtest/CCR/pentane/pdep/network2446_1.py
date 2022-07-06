species(
    label = 'C[CH][CH]C[CH][C]=O(4251)',
    structure = SMILES('C[CH][CH]C[CH][C]=O'),
    E0 = (441.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1855,455,950,458.984,2997.16,4000],'cm^-1')),
        HinderedRotor(inertia=(0.037149,'amu*angstrom^2'), symmetry=1, barrier=(4.85483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037149,'amu*angstrom^2'), symmetry=1, barrier=(4.85483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037149,'amu*angstrom^2'), symmetry=1, barrier=(4.85483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037149,'amu*angstrom^2'), symmetry=1, barrier=(4.85483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037149,'amu*angstrom^2'), symmetry=1, barrier=(4.85483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36113,0.066249,-9.0788e-05,8.42002e-08,-3.10858e-11,53215,31.4754], Tmin=(100,'K'), Tmax=(858.601,'K')), NASAPolynomial(coeffs=[2.08663,0.0441922,-1.96253e-05,3.61042e-09,-2.43553e-13,53778.9,32.0948], Tmin=(858.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJCHO) + radical(CCCJ=O)"""),
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
    label = 'C[CH][CH]CC=C=O(4248)',
    structure = SMILES('C[CH][CH]CC=C=O'),
    E0 = (239.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,1534.87,1534.88],'cm^-1')),
        HinderedRotor(inertia=(0.197241,'amu*angstrom^2'), symmetry=1, barrier=(4.53496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197297,'amu*angstrom^2'), symmetry=1, barrier=(4.53625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00271424,'amu*angstrom^2'), symmetry=1, barrier=(4.53718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197318,'amu*angstrom^2'), symmetry=1, barrier=(4.53673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46779,0.0419894,-1.81794e-05,3.05303e-09,-1.68411e-13,28810.7,22.3028], Tmin=(100,'K'), Tmax=(2999.39,'K')), NASAPolynomial(coeffs=[48.2198,-0.00812339,2.09681e-06,-3.90076e-10,2.99185e-14,-2339.22,-249.931], Tmin=(2999.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH][C]=O(9632)',
    structure = SMILES('[CH2][CH][C]=O'),
    E0 = (334.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.0670133,'amu*angstrom^2'), symmetry=1, barrier=(7.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00540518,'amu*angstrom^2'), symmetry=1, barrier=(7.51488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.78792,0.0268854,-2.93724e-05,1.83573e-08,-4.61895e-12,40303.1,15.8795], Tmin=(100,'K'), Tmax=(968.231,'K')), NASAPolynomial(coeffs=[6.80976,0.01027,-3.63122e-06,6.33225e-10,-4.24838e-14,39524.3,-3.39371], Tmin=(968.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJCHO) + radical(CJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][CH]C(3874)',
    structure = SMILES('[CH][CH]C'),
    E0 = (522.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1852.85,1853.13,1853.16],'cm^-1')),
        HinderedRotor(inertia=(0.0369856,'amu*angstrom^2'), symmetry=1, barrier=(8.28142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370949,'amu*angstrom^2'), symmetry=1, barrier=(8.27684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19279,0.0163284,-1.76328e-06,-3.55309e-09,1.20622e-12,62877.3,14.3587], Tmin=(100,'K'), Tmax=(1426.65,'K')), NASAPolynomial(coeffs=[5.59458,0.0149965,-6.04303e-06,1.1011e-09,-7.44871e-14,61642.3,-0.00873183], Tmin=(1426.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(CCJ2_triplet)"""),
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
    label = '[CH][CH]C[CH][C]=O(9824)',
    structure = SMILES('[CH][CH]C[CH][C]=O'),
    E0 = (719.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,1855,455,950,214.637,906.452,1176.38,1616.17,1985.53],'cm^-1')),
        HinderedRotor(inertia=(0.104485,'amu*angstrom^2'), symmetry=1, barrier=(3.0842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104485,'amu*angstrom^2'), symmetry=1, barrier=(3.0842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104485,'amu*angstrom^2'), symmetry=1, barrier=(3.0842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104485,'amu*angstrom^2'), symmetry=1, barrier=(3.0842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72568,0.0537284,-7.12007e-05,5.75478e-08,-1.8864e-11,86586.3,26.3634], Tmin=(100,'K'), Tmax=(859.917,'K')), NASAPolynomial(coeffs=[6.33535,0.0257354,-1.09443e-05,1.97419e-09,-1.31908e-13,86035.7,6.22837], Tmin=(859.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C[C][CH]C[CH][C]=O(10313)',
    structure = SMILES('C[C][CH]C[CH][C]=O'),
    E0 = (695.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09384,0.0692637,-9.13036e-05,7.5639e-08,-2.57198e-11,83748.8,29.6061], Tmin=(100,'K'), Tmax=(834.003,'K')), NASAPolynomial(coeffs=[6.22772,0.0359921,-1.59075e-05,2.93644e-09,-1.9931e-13,83193.2,7.57323], Tmin=(834.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][C]C[CH][C]=O(10314)',
    structure = SMILES('C[CH][C]C[CH][C]=O'),
    E0 = (695.488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06415,0.0686879,-8.57072e-05,6.66422e-08,-2.15731e-11,83749.6,29.8852], Tmin=(100,'K'), Tmax=(821.286,'K')), NASAPolynomial(coeffs=[7.2909,0.0340934,-1.47292e-05,2.69988e-09,-1.83016e-13,82870.8,1.94695], Tmin=(821.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH]C[C][C]=O(10315)',
    structure = SMILES('C[CH][CH]C[C][C]=O'),
    E0 = (722.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,247.496,859.297,1246.71,1853],'cm^-1')),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14077,0.0757163,-0.000127953,1.26598e-07,-4.71935e-11,86977.3,32.4101], Tmin=(100,'K'), Tmax=(870.25,'K')), NASAPolynomial(coeffs=[1.2671,0.0443492,-2.08223e-05,3.87836e-09,-2.61356e-13,88121.1,38.5161], Tmin=(870.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C1CC1[C]=O(10303)',
    structure = SMILES('C[CH]C1CC1[C]=O'),
    E0 = (226.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55418,0.0424575,5.55433e-06,-3.30429e-08,1.427e-11,27322.8,27.0614], Tmin=(100,'K'), Tmax=(1029.51,'K')), NASAPolynomial(coeffs=[11.76,0.0284561,-1.14202e-05,2.15116e-09,-1.53308e-13,23862,-29.0751], Tmin=(1029.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Cs_S) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C[CH]CC=C[C]=O(10316)',
    structure = SMILES('C[CH]CC=C[C]=O'),
    E0 = (191.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00827,0.0721189,-9.24406e-05,7.93271e-08,-2.93044e-11,23107.1,28.399], Tmin=(100,'K'), Tmax=(743.49,'K')), NASAPolynomial(coeffs=[5.34052,0.0427122,-2.08074e-05,4.06217e-09,-2.86407e-13,22631.5,9.91618], Tmin=(743.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(RCCJC) + radical(C=CCJ=O)"""),
)

species(
    label = 'C[CH]C=CC[C]=O(4247)',
    structure = SMILES('C[CH]C=CC[C]=O'),
    E0 = (149.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,180,1087.05],'cm^-1')),
        HinderedRotor(inertia=(0.0187901,'amu*angstrom^2'), symmetry=1, barrier=(15.7432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684932,'amu*angstrom^2'), symmetry=1, barrier=(15.7479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188262,'amu*angstrom^2'), symmetry=1, barrier=(15.7399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122308,'amu*angstrom^2'), symmetry=1, barrier=(15.7516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36151,0.0467316,-3.90031e-06,-2.24189e-08,9.72247e-12,18056.7,25.5773], Tmin=(100,'K'), Tmax=(1120.66,'K')), NASAPolynomial(coeffs=[12.9783,0.0295868,-1.35031e-05,2.65778e-09,-1.91472e-13,13925.8,-38.6038], Tmin=(1120.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C[CH][CH][C]=O(10022)',
    structure = SMILES('C[CH]C=C[CH][C]=O'),
    E0 = (311.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,570.482,570.653],'cm^-1')),
        HinderedRotor(inertia=(0.106538,'amu*angstrom^2'), symmetry=1, barrier=(24.5795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000517607,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106493,'amu*angstrom^2'), symmetry=1, barrier=(24.5681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106445,'amu*angstrom^2'), symmetry=1, barrier=(24.5784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48667,0.039376,2.1295e-05,-5.40161e-08,2.21245e-11,37556.2,25.5063], Tmin=(100,'K'), Tmax=(1028.65,'K')), NASAPolynomial(coeffs=[15.2159,0.0235917,-1.05212e-05,2.14113e-09,-1.60646e-13,32742.2,-50.787], Tmin=(1028.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_S) + radical(CCJC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=CC[CH][C]=O(5329)',
    structure = SMILES('[CH2]C=CC[CH][C]=O'),
    E0 = (319.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.778956,'amu*angstrom^2'), symmetry=1, barrier=(17.9097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403324,'amu*angstrom^2'), symmetry=1, barrier=(50.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779054,'amu*angstrom^2'), symmetry=1, barrier=(17.912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968321,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3657.02,'J/mol'), sigma=(6.16401,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.22 K, Pc=35.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907921,0.0595326,-4.70532e-05,1.94918e-08,-3.25348e-12,38588.9,28.1278], Tmin=(100,'K'), Tmax=(1429.35,'K')), NASAPolynomial(coeffs=[14.097,0.0226237,-8.32024e-06,1.42646e-09,-9.37883e-14,34818.5,-40.2132], Tmin=(1429.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'CH2CHCO(3668)',
    structure = SMILES('C=C[C]=O'),
    E0 = (83.3963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.27992,'amu*angstrom^2'), symmetry=1, barrier=(29.4278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08334,0.0153205,6.54259e-06,-1.77535e-08,7.39369e-12,10067.5,11.8951], Tmin=(100,'K'), Tmax=(1002.99,'K')), NASAPolynomial(coeffs=[6.9784,0.0111885,-4.32948e-06,8.06749e-10,-5.75264e-14,8712.67,-9.76694], Tmin=(1002.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH][CH][CH][C]=O(10317)',
    structure = SMILES('C[CH][CH][CH]C=[C][O]'),
    E0 = (638.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3010,987.5,1337.5,450,1655,1685,370,180,1115.72,1115.78,1115.85],'cm^-1')),
        HinderedRotor(inertia=(0.181327,'amu*angstrom^2'), symmetry=1, barrier=(4.16907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00472112,'amu*angstrom^2'), symmetry=1, barrier=(4.17203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18076,'amu*angstrom^2'), symmetry=1, barrier=(4.15602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00473836,'amu*angstrom^2'), symmetry=1, barrier=(4.16996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69464,0.0546385,-5.51821e-05,4.2767e-08,-1.53176e-11,76854.9,31.4366], Tmin=(100,'K'), Tmax=(764.27,'K')), NASAPolynomial(coeffs=[3.86453,0.0393853,-1.75976e-05,3.31129e-09,-2.2914e-13,76637.1,22.2959], Tmin=(764.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(RCCJCC) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][C]=O(9128)',
    structure = SMILES('[CH2][CH][CH]C[CH][C]=O'),
    E0 = (646.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1855,455,950,401.252,3128.77,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38728,0.0671154,-0.000101705,9.70352e-08,-3.57488e-11,77898.2,33.1678], Tmin=(100,'K'), Tmax=(870.5,'K')), NASAPolynomial(coeffs=[1.93353,0.0418045,-1.88017e-05,3.45454e-09,-2.31828e-13,78667,35.5702], Tmin=(870.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C[CH][CH][C]=O(10318)',
    structure = SMILES('C[CH]C[CH]C=[C][O]'),
    E0 = (443.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,434.359,434.36,434.362,1784.18],'cm^-1')),
        HinderedRotor(inertia=(0.0634591,'amu*angstrom^2'), symmetry=1, barrier=(8.496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0017424,'amu*angstrom^2'), symmetry=1, barrier=(19.7833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0423545,'amu*angstrom^2'), symmetry=1, barrier=(95.6752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311331,'amu*angstrom^2'), symmetry=1, barrier=(41.6819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29796,0.0558319,-3.6704e-05,1.25855e-08,-1.80123e-12,53487.6,30.2942], Tmin=(100,'K'), Tmax=(1573.21,'K')), NASAPolynomial(coeffs=[11.4476,0.0300257,-1.20988e-05,2.15876e-09,-1.44303e-13,50294.1,-23.271], Tmin=(1573.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][CH][CH]C[C]=O(4250)',
    structure = SMILES('C[CH][CH][CH]C[C]=O'),
    E0 = (474.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1855,455,950,209.497,2798.29,3080.95],'cm^-1')),
        HinderedRotor(inertia=(0.0185402,'amu*angstrom^2'), symmetry=1, barrier=(3.01885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0185402,'amu*angstrom^2'), symmetry=1, barrier=(3.01885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0185402,'amu*angstrom^2'), symmetry=1, barrier=(3.01885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0185402,'amu*angstrom^2'), symmetry=1, barrier=(3.01885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0185402,'amu*angstrom^2'), symmetry=1, barrier=(3.01885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39219,0.0673961,-9.96537e-05,9.72889e-08,-3.69555e-11,57105.7,34.0987], Tmin=(100,'K'), Tmax=(855.507,'K')), NASAPolynomial(coeffs=[0.876135,0.0462345,-2.1216e-05,3.95483e-09,-2.6839e-13,58056.7,41.5499], Tmin=(855.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC) + radical(RCCJCC) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C[CH]C[CH][C]=O(5322)',
    structure = SMILES('[CH2]C[CH]C[CH][C]=O'),
    E0 = (452.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,191.599,1487.3,4000],'cm^-1')),
        HinderedRotor(inertia=(0.12169,'amu*angstrom^2'), symmetry=1, barrier=(3.09029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12169,'amu*angstrom^2'), symmetry=1, barrier=(3.09029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12169,'amu*angstrom^2'), symmetry=1, barrier=(3.09029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12169,'amu*angstrom^2'), symmetry=1, barrier=(3.09029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12169,'amu*angstrom^2'), symmetry=1, barrier=(3.09029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20298,0.0667101,-8.144e-05,6.72404e-08,-2.34182e-11,54522.6,31.0957], Tmin=(100,'K'), Tmax=(819.487,'K')), NASAPolynomial(coeffs=[5.02648,0.0399478,-1.76288e-05,3.26841e-09,-2.23009e-13,54167.9,15.07], Tmin=(819.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[CH][CH][CH][C]=O(9754)',
    structure = SMILES('CC[CH][CH]C=[C][O]'),
    E0 = (443.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,180,724.128,724.131,724.151],'cm^-1')),
        HinderedRotor(inertia=(0.00516824,'amu*angstrom^2'), symmetry=1, barrier=(1.92307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836222,'amu*angstrom^2'), symmetry=1, barrier=(1.92264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836423,'amu*angstrom^2'), symmetry=1, barrier=(1.9231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836431,'amu*angstrom^2'), symmetry=1, barrier=(1.92312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54721,0.0539207,-3.3969e-05,1.11053e-08,-1.5421e-12,53476.9,29.2191], Tmin=(100,'K'), Tmax=(1566.39,'K')), NASAPolynomial(coeffs=[9.81206,0.0328153,-1.37582e-05,2.50346e-09,-1.69231e-13,50887.7,-14.3631], Tmin=(1566.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]CC[CH][C]=O(5321)',
    structure = SMILES('[CH2][CH]CC[CH][C]=O'),
    E0 = (452.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,221.051,421.53,3442.46],'cm^-1')),
        HinderedRotor(inertia=(0.0260358,'amu*angstrom^2'), symmetry=1, barrier=(0.788456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260358,'amu*angstrom^2'), symmetry=1, barrier=(0.788456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260358,'amu*angstrom^2'), symmetry=1, barrier=(0.788456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260358,'amu*angstrom^2'), symmetry=1, barrier=(0.788456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260358,'amu*angstrom^2'), symmetry=1, barrier=(0.788456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18132,0.0660295,-7.54261e-05,5.76226e-08,-1.89677e-11,54523.1,31.3467], Tmin=(100,'K'), Tmax=(800.678,'K')), NASAPolynomial(coeffs=[6.06893,0.0380867,-1.6473e-05,3.03737e-09,-2.07184e-13,53853.4,9.55902], Tmin=(800.678,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH][CH][CH]C=O(10319)',
    structure = SMILES('C[CH][CH][CH]C=C[O]'),
    E0 = (398.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29356,0.0550838,-3.54871e-05,1.19991e-08,-1.69484e-12,48041.9,29.7476], Tmin=(100,'K'), Tmax=(1594.53,'K')), NASAPolynomial(coeffs=[11.3689,0.0298091,-1.17108e-05,2.05835e-09,-1.36262e-13,44828.9,-23.5608], Tmin=(1594.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(RCCJCC) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][CH][CH]CC[C]=O(4252)',
    structure = SMILES('[CH2][CH][CH]CC[C]=O'),
    E0 = (479.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,270.054,1764.02,4000],'cm^-1')),
        HinderedRotor(inertia=(0.020017,'amu*angstrom^2'), symmetry=1, barrier=(3.2067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.020017,'amu*angstrom^2'), symmetry=1, barrier=(3.2067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.020017,'amu*angstrom^2'), symmetry=1, barrier=(3.2067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.020017,'amu*angstrom^2'), symmetry=1, barrier=(3.2067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.020017,'amu*angstrom^2'), symmetry=1, barrier=(3.2067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3572.62,'J/mol'), sigma=(6.29932,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.04 K, Pc=32.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23894,0.0732996,-0.000118599,1.18895e-07,-4.51957e-11,57751.6,33.9385], Tmin=(100,'K'), Tmax=(864.247,'K')), NASAPolynomial(coeffs=[0.110065,0.0482266,-2.2497e-05,4.19909e-09,-2.84107e-13,59078.2,45.7662], Tmin=(864.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]C=O(5551)',
    structure = SMILES('[CH2][CH][CH]CC=C[O]'),
    E0 = (462.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,991.213,991.493],'cm^-1')),
        HinderedRotor(inertia=(0.162799,'amu*angstrom^2'), symmetry=1, barrier=(3.74306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162033,'amu*angstrom^2'), symmetry=1, barrier=(3.72546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00536787,'amu*angstrom^2'), symmetry=1, barrier=(3.73664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00535731,'amu*angstrom^2'), symmetry=1, barrier=(3.73103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41551,0.0577798,-4.39124e-05,1.92192e-08,-3.66946e-12,55746.8,32.1999], Tmin=(100,'K'), Tmax=(1188.33,'K')), NASAPolynomial(coeffs=[8.37237,0.0343626,-1.43534e-05,2.63623e-09,-1.80753e-13,54093.4,-2.56335], Tmin=(1188.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2049,0.0625364,-5.72642e-05,3.16551e-08,-7.47846e-12,53832.1,30.8], Tmin=(100,'K'), Tmax=(999.296,'K')), NASAPolynomial(coeffs=[8.4989,0.0333397,-1.3438e-05,2.41688e-09,-1.63698e-13,52374.4,-4.38418], Tmin=(999.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C1CC=C1[O](10320)',
    structure = SMILES('C[CH]C1C[CH]C1=O'),
    E0 = (192.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90348,0.0260053,6.68967e-05,-1.04047e-07,4.10489e-11,23279.7,25.5224], Tmin=(100,'K'), Tmax=(967.705,'K')), NASAPolynomial(coeffs=[14.4493,0.0244225,-8.57938e-06,1.63671e-09,-1.2331e-13,18497.6,-46.7548], Tmin=(967.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJCC=O) + radical(CCJC=O)"""),
)

species(
    label = 'C[CH]C=CC=C[O](10321)',
    structure = SMILES('CC=C[CH]C=C[O]'),
    E0 = (87.1191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27304,0.0401786,3.92183e-05,-8.6927e-08,3.83781e-11,10594.3,25.3567], Tmin=(100,'K'), Tmax=(939.29,'K')), NASAPolynomial(coeffs=[18.38,0.0173919,-4.34176e-06,7.34531e-10,-5.68058e-14,5172.19,-67.8589], Tmin=(939.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.1191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC=CC[C]=[C][O](10038)',
    structure = SMILES('CC=CC[C][C]=O'),
    E0 = (449.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,318.548,318.548,318.549],'cm^-1')),
        HinderedRotor(inertia=(0.14232,'amu*angstrom^2'), symmetry=1, barrier=(10.2481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142319,'amu*angstrom^2'), symmetry=1, barrier=(10.2481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142319,'amu*angstrom^2'), symmetry=1, barrier=(10.248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142319,'amu*angstrom^2'), symmetry=1, barrier=(10.2481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23752,0.0649474,-6.96586e-05,4.5146e-08,-1.25269e-11,54103.9,26.8144], Tmin=(100,'K'), Tmax=(856.174,'K')), NASAPolynomial(coeffs=[7.85775,0.0340181,-1.54715e-05,2.95305e-09,-2.06772e-13,52970.3,-4.09639], Tmin=(856.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH]C[C]=C[O](10322)',
    structure = SMILES('C[CH][CH]C[C]=C[O]'),
    E0 = (495.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,401.184,401.184,2151.71,2151.72],'cm^-1')),
        HinderedRotor(inertia=(0.00357064,'amu*angstrom^2'), symmetry=1, barrier=(11.7312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102713,'amu*angstrom^2'), symmetry=1, barrier=(11.7312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60963,'amu*angstrom^2'), symmetry=1, barrier=(69.6274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102714,'amu*angstrom^2'), symmetry=1, barrier=(11.7312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36187,0.0613006,-5.64328e-05,3.32187e-08,-8.70443e-12,59667.1,31.0906], Tmin=(100,'K'), Tmax=(888.86,'K')), NASAPolynomial(coeffs=[6.66547,0.0374337,-1.61561e-05,3.01027e-09,-2.08058e-13,58724.2,6.12864], Tmin=(888.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]CC[C]=[C][O](10323)',
    structure = SMILES('C[CH]CC[C][C]=O'),
    E0 = (527.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901789,0.0750558,-0.000103339,8.96183e-08,-3.15798e-11,63603.6,30.7052], Tmin=(100,'K'), Tmax=(831.593,'K')), NASAPolynomial(coeffs=[5.49172,0.0404713,-1.83978e-05,3.43792e-09,-2.34741e-13,63272.6,12.008], Tmin=(831.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH]C[C]=[C]O(10324)',
    structure = SMILES('C[CH][CH]C[C]=[C]O'),
    E0 = (593.616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1670,1700,300,440,208.28,808.28,1622.47],'cm^-1')),
        HinderedRotor(inertia=(0.141986,'amu*angstrom^2'), symmetry=1, barrier=(3.54042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141986,'amu*angstrom^2'), symmetry=1, barrier=(3.54042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141986,'amu*angstrom^2'), symmetry=1, barrier=(3.54042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141986,'amu*angstrom^2'), symmetry=1, barrier=(3.54042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141986,'amu*angstrom^2'), symmetry=1, barrier=(3.54042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950697,0.0726612,-9.67094e-05,7.99232e-08,-2.67683e-11,71499.9,34.274], Tmin=(100,'K'), Tmax=(858.183,'K')), NASAPolynomial(coeffs=[6.41736,0.0367533,-1.57204e-05,2.84908e-09,-1.90915e-13,70945.6,10.9738], Tmin=(858.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJCC) + radical(RCCJC) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'CC[CH]C[C][C]=O(9751)',
    structure = SMILES('CC[CH]C[C][C]=O'),
    E0 = (527.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.937833,0.0755488,-0.000108606,9.81234e-08,-3.54841e-11,63602.5,30.4037], Tmin=(100,'K'), Tmax=(837.553,'K')), NASAPolynomial(coeffs=[4.41443,0.0423955,-1.95915e-05,3.67821e-09,-2.51351e-13,63600.6,17.7128], Tmin=(837.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH][CH]C=[C]O(10325)',
    structure = SMILES('C[CH]C=C[CH][C]O'),
    E0 = (483.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.305873,0.0730937,-6.412e-05,2.83446e-08,-4.97795e-12,58351.2,25.5905], Tmin=(100,'K'), Tmax=(1369.4,'K')), NASAPolynomial(coeffs=[17.3954,0.0231753,-9.44075e-06,1.72492e-09,-1.18211e-13,53670.7,-62.229], Tmin=(1369.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=CCJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2][CH][CH]CC=[C]O(10326)',
    structure = SMILES('[CH2][CH][CH]CC=[C]O'),
    E0 = (561.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,298.584,455.383,3148.11],'cm^-1')),
        HinderedRotor(inertia=(0.0750225,'amu*angstrom^2'), symmetry=1, barrier=(1.91776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750225,'amu*angstrom^2'), symmetry=1, barrier=(1.91776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750225,'amu*angstrom^2'), symmetry=1, barrier=(1.91776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750225,'amu*angstrom^2'), symmetry=1, barrier=(1.91776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750225,'amu*angstrom^2'), symmetry=1, barrier=(1.91776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05594,0.0683475,-8.06386e-05,6.04012e-08,-1.90523e-11,67577.7,35.2113], Tmin=(100,'K'), Tmax=(820.465,'K')), NASAPolynomial(coeffs=[7.27725,0.0351771,-1.48037e-05,2.68867e-09,-1.81606e-13,66652.5,7.01068], Tmin=(820.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
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
    E0 = (441.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (441.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (912.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (930.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (911.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (907.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (907.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (934.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (444.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (464.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (464.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (530.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (539.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (546.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (605.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (850.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (858.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (602.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (632.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (594.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (595.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (596.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (594.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (553.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (536.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (604.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (450.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (505.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (675.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (637.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (699.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (669.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (743.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (642.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (567.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (657.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['HCCO(2227)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH][CH]CC=C=O(4248)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH][C]=O(9632)', '[CH][CH]C(3874)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH][CH]C(3856)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(17)', '[CH][CH]C[CH][C]=O(9824)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C[C][CH]C[CH][C]=O(10313)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C[CH][C]C[CH][C]=O(10314)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C[CH][CH]C[C][C]=O(10315)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH]C1CC1[C]=O(10303)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH]CC=C[C]=O(10316)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH]C=CC[C]=O(4247)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', 'CC=C[CH][CH][C]=O(10022)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C][O](6861)', 'm1_allyl(186)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2CHCO(3668)', '[CH][CH]C(3874)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C[CH][CH][CH][CH][C]=O(10317)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH]C[CH][CH][C]=O(10318)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[CH][CH][CH]C[C]=O(4250)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C[CH]C[CH][C]=O(5322)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC[CH][CH][CH][C]=O(9754)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.86945e+08,'s^-1'), n=1.27834, Ea=(151.688,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC[CH][C]=O(5321)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH][CH][CH][CH]C=O(10319)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.21412e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH][CH]CC[C]=O(4252)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH]C[CH]C=O(5551)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6Hall;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH]C)[CH][C]=O(4240)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH]C1CC=C1[O](10320)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[CH][CH]C[CH][C]=O(4251)'],
    products = ['C[CH]C=CC=C[O](10321)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', 'CC=CC[C]=[C][O](10038)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['HCCO(2227)', '[CH2][CH][CH]C(3856)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH][CH]C[C]=C[O](10322)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[CH]CC[C]=[C][O](10323)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH][CH]C[C]=[C]O(10324)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CC[CH]C[C][C]=O(9751)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(893531,'s^-1'), n=1.72475, Ea=(114.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;Cs_H_out_H/NonDeC] + [R4Hall;Cd_rad_out;Cs_H_out_1H] for rate rule [R4HJ_2;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH][CH][CH]C=[C]O(10325)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH][CH]CC=[C]O(10326)'],
    products = ['C[CH][CH]C[CH][C]=O(4251)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(30.079,'s^-1'), n=2.77074, Ea=(96.0343,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;O_H_out] + [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7Hall;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2446',
    isomers = [
        'C[CH][CH]C[CH][C]=O(4251)',
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
    label = 'PDepNetwork #2446',
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


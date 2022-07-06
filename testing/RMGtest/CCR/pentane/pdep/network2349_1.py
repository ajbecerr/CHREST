species(
    label = '[CH2]C([CH][C]=O)CO[O](2913)',
    structure = SMILES('[CH2]C([CH][C]=O)CO[O]'),
    E0 = (336.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,2421.64,2422.4],'cm^-1')),
        HinderedRotor(inertia=(1.30385,'amu*angstrom^2'), symmetry=1, barrier=(69.3665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153415,'amu*angstrom^2'), symmetry=1, barrier=(8.16968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154129,'amu*angstrom^2'), symmetry=1, barrier=(8.1748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153437,'amu*angstrom^2'), symmetry=1, barrier=(8.17244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30711,'amu*angstrom^2'), symmetry=1, barrier=(69.3617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.423198,0.0832783,-0.000117753,9.25598e-08,-2.86563e-11,40575.5,33.0248], Tmin=(100,'K'), Tmax=(896.328,'K')), NASAPolynomial(coeffs=[10.3546,0.0294991,-1.19253e-05,2.07417e-09,-1.34732e-13,39175.1,-11.682], Tmin=(896.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH]C([C]=O)CO[O](10551)',
    structure = SMILES('[CH2][CH]C([C]=O)CO[O]'),
    E0 = (370.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,286.472,1412.68],'cm^-1')),
        HinderedRotor(inertia=(0.00205412,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916748,'amu*angstrom^2'), symmetry=1, barrier=(5.33883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.091676,'amu*angstrom^2'), symmetry=1, barrier=(5.33887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916751,'amu*angstrom^2'), symmetry=1, barrier=(5.33884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0916764,'amu*angstrom^2'), symmetry=1, barrier=(5.3389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05703,0.0729496,-7.99084e-05,3.24261e-08,6.83672e-12,44704.1,34.0435], Tmin=(100,'K'), Tmax=(571.934,'K')), NASAPolynomial(coeffs=[8.83331,0.0334773,-1.54987e-05,2.93998e-09,-2.03477e-13,43570.6,-1.26008], Tmin=(571.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCC=O) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[O]OCC[CH][CH][C]=O(10552)',
    structure = SMILES('[O][C]=C[CH]CCO[O]'),
    E0 = (336.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,341.333,698.966,698.97,2337.41],'cm^-1')),
        HinderedRotor(inertia=(0.532714,'amu*angstrom^2'), symmetry=1, barrier=(12.2481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.532688,'amu*angstrom^2'), symmetry=1, barrier=(12.2476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.532701,'amu*angstrom^2'), symmetry=1, barrier=(12.2478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208437,'amu*angstrom^2'), symmetry=1, barrier=(72.2586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894674,0.0704696,-7.23641e-05,4.0995e-08,-9.59119e-12,40630.6,30.9928], Tmin=(100,'K'), Tmax=(1020.92,'K')), NASAPolynomial(coeffs=[11.1181,0.0304125,-1.35077e-05,2.56006e-09,-1.79034e-13,38543.2,-18.5406], Tmin=(1020.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
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
    label = '[CH2]C([CH][C]=O)C[O](10553)',
    structure = SMILES('[CH2]C([CH][C]=O)C[O]'),
    E0 = (338.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,1262.83,1262.83],'cm^-1')),
        HinderedRotor(inertia=(0.0860189,'amu*angstrom^2'), symmetry=1, barrier=(1.97774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860206,'amu*angstrom^2'), symmetry=1, barrier=(1.97778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860188,'amu*angstrom^2'), symmetry=1, barrier=(1.97774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860025,'amu*angstrom^2'), symmetry=1, barrier=(1.97737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2173,0.0669511,-9.32998e-05,7.8748e-08,-2.64624e-11,40809.8,29.5147], Tmin=(100,'K'), Tmax=(875.073,'K')), NASAPolynomial(coeffs=[6.04537,0.0327786,-1.39768e-05,2.5109e-09,-1.66735e-13,40428.3,9.51416], Tmin=(875.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[O]OC[CH][CH][C]=O(10554)',
    structure = SMILES('[O][C]=C[CH]CO[O]'),
    E0 = (336.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,292.73,293.038,293.617],'cm^-1')),
        HinderedRotor(inertia=(0.270186,'amu*angstrom^2'), symmetry=1, barrier=(16.5445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.50758,'amu*angstrom^2'), symmetry=1, barrier=(31.0115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2187,'amu*angstrom^2'), symmetry=1, barrier=(74.8171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09561,0.0652133,-7.64351e-05,4.68025e-08,-1.14485e-11,40574.1,26.3311], Tmin=(100,'K'), Tmax=(994.434,'K')), NASAPolynomial(coeffs=[12.1714,0.0206619,-9.23356e-06,1.75039e-09,-1.22404e-13,38371.3,-27.0412], Tmin=(994.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ) + radical(C=CCJCO) + radical(C=CJO)"""),
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
    label = '[CH2]C([C][C]=O)CO[O](10555)',
    structure = SMILES('[CH2]C([C][C]=O)CO[O]'),
    E0 = (617.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.21555,0.0925874,-0.000154328,1.34144e-07,-4.43991e-11,74337.2,33.9145], Tmin=(100,'K'), Tmax=(895.617,'K')), NASAPolynomial(coeffs=[9.48987,0.0297356,-1.31693e-05,2.35338e-09,-1.53482e-13,73535.5,-5.00784], Tmin=(895.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C([CH][C]=O)CO[O](10556)',
    structure = SMILES('[CH]C([CH][C]=O)CO[O]'),
    E0 = (579.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.395914,0.0841411,-0.000123045,9.76889e-08,-3.06834e-11,69818.3,32.3456], Tmin=(100,'K'), Tmax=(853.5,'K')), NASAPolynomial(coeffs=[11.166,0.0272982,-1.19541e-05,2.17432e-09,-1.45552e-13,68211.8,-16.5491], Tmin=(853.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OCC1CC1[C]=O(10557)',
    structure = SMILES('[O]OCC1CC1[C]=O'),
    E0 = (119.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844017,0.0622175,-5.07956e-05,2.13771e-08,-3.61685e-12,14465.8,29.0382], Tmin=(100,'K'), Tmax=(1408.35,'K')), NASAPolynomial(coeffs=[14.6438,0.0230233,-9.0505e-06,1.61621e-09,-1.09017e-13,10578.9,-42.2628], Tmin=(1408.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(ROOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1COOC1[C]=O(10558)',
    structure = SMILES('[CH2]C1COOC1[C]=O'),
    E0 = (91.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32043,0.0403085,3.86561e-05,-9.57207e-08,4.62046e-11,11079.1,28.1415], Tmin=(100,'K'), Tmax=(878.615,'K')), NASAPolynomial(coeffs=[19.4256,0.00978692,2.15077e-06,-7.84778e-10,6.02286e-14,5894.18,-68.2635], Tmin=(878.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C1COOC1(10559)',
    structure = SMILES('O=[C][CH]C1COOC1'),
    E0 = (63.6682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2356,0.0485589,5.14901e-06,-5.29616e-08,2.88881e-11,7768.78,25.0075], Tmin=(100,'K'), Tmax=(869.769,'K')), NASAPolynomial(coeffs=[15.9625,0.0167663,-1.99527e-06,1.64384e-11,6.90318e-15,3847.73,-51.8003], Tmin=(869.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.6682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1COC1[C]=O(10560)',
    structure = SMILES('[CH2]C1COC1[C]=O'),
    E0 = (119.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90481,0.0305747,4.38194e-05,-9.14384e-08,4.29463e-11,14434.9,25.2589], Tmin=(100,'K'), Tmax=(867.736,'K')), NASAPolynomial(coeffs=[15.5262,0.0115086,1.19411e-06,-6.20403e-10,5.08483e-14,10424.8,-48.0094], Tmin=(867.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C1COC1(10561)',
    structure = SMILES('O=[C][CH]C1COC1'),
    E0 = (91.7692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82852,0.0387077,1.08163e-05,-4.95011e-08,2.60743e-11,11124.3,21.4021], Tmin=(100,'K'), Tmax=(848.969,'K')), NASAPolynomial(coeffs=[12.0677,0.0184815,-2.9488e-06,1.80183e-10,-2.43151e-15,8376.06,-32.2661], Tmin=(848.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.7692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CCJCHO) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([CH2])=C[C]=O(9813)',
    structure = SMILES('[CH2]C([CH2])=C[C]=O'),
    E0 = (307.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.12597,'amu*angstrom^2'), symmetry=1, barrier=(25.8882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12599,'amu*angstrom^2'), symmetry=1, barrier=(25.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.126,'amu*angstrom^2'), symmetry=1, barrier=(25.889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50488,0.0572345,-5.84403e-05,3.09075e-08,-6.68244e-12,37035.6,20.0483], Tmin=(100,'K'), Tmax=(1099.74,'K')), NASAPolynomial(coeffs=[11.0466,0.022529,-1.11032e-05,2.21137e-09,-1.58996e-13,34936.9,-26.8919], Tmin=(1099.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=C(C[C]=O)CO[O](2910)',
    structure = SMILES('C=C(C[C]=O)CO[O]'),
    E0 = (98.1407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,673.452,4000],'cm^-1')),
        HinderedRotor(inertia=(0.392217,'amu*angstrom^2'), symmetry=1, barrier=(9.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799782,'amu*angstrom^2'), symmetry=1, barrier=(18.3886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.80006,'amu*angstrom^2'), symmetry=1, barrier=(18.395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799839,'amu*angstrom^2'), symmetry=1, barrier=(18.3899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36676,0.064003,-5.77296e-05,2.83881e-08,-6.03979e-12,11893.3,27.981], Tmin=(100,'K'), Tmax=(1073.47,'K')), NASAPolynomial(coeffs=[9.04767,0.0353826,-1.77377e-05,3.55198e-09,-2.55783e-13,10244.3,-9.61956], Tmin=(1073.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.1407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC(=C[C]=O)CO[O](10562)',
    structure = SMILES('CC(=C[C]=O)CO[O]'),
    E0 = (74.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.221715,0.0941092,-0.00015726,1.45307e-07,-5.2376e-11,9094.83,29.358], Tmin=(100,'K'), Tmax=(826.033,'K')), NASAPolynomial(coeffs=[7.12861,0.0396581,-2.02386e-05,3.93688e-09,-2.73332e-13,8670.38,1.69397], Tmin=(826.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C(=C[C]=O)COO(10563)',
    structure = SMILES('[CH2]C(=C[C]=O)COO'),
    E0 = (74.0705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118715,0.0931091,-0.000135702,1.10462e-07,-3.67879e-11,9041.14,29.2619], Tmin=(100,'K'), Tmax=(730.704,'K')), NASAPolynomial(coeffs=[10.5229,0.0361443,-1.8743e-05,3.73366e-09,-2.65601e-13,7520.93,-17.6662], Tmin=(730.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.0705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C(=C[C]=O)CO[O](10564)',
    structure = SMILES('[CH2]C(=C[C]=O)CO[O]'),
    E0 = (226.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.574214,'amu*angstrom^2'), symmetry=1, barrier=(13.2023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909566,'amu*angstrom^2'), symmetry=1, barrier=(20.9127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.573911,'amu*angstrom^2'), symmetry=1, barrier=(13.1953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574217,'amu*angstrom^2'), symmetry=1, barrier=(13.2024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.236997,0.0914031,-0.000147805,1.30497e-07,-4.5607e-11,27317.7,29.5852], Tmin=(100,'K'), Tmax=(806.221,'K')), NASAPolynomial(coeffs=[9.34289,0.0337435,-1.73051e-05,3.38384e-09,-2.36199e-13,26255.1,-9.86843], Tmin=(806.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]OCC=C[C]=O(10565)',
    structure = SMILES('[O]OC[CH]C=C=O'),
    E0 = (78.9155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.79974,'amu*angstrom^2'), symmetry=1, barrier=(41.3797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.57994,'amu*angstrom^2'), symmetry=1, barrier=(13.334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.54626,'amu*angstrom^2'), symmetry=1, barrier=(81.5355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13866,0.0674966,-9.2016e-05,7.01733e-08,-2.18129e-11,9590.17,22.5155], Tmin=(100,'K'), Tmax=(782.865,'K')), NASAPolynomial(coeffs=[9.26787,0.0259592,-1.24257e-05,2.3936e-09,-1.67288e-13,8317.4,-14.713], Tmin=(782.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.9155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=CCJCO)"""),
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
    label = 'C=CC=[C][O](9589)',
    structure = SMILES('[CH2]C=C[C]=O'),
    E0 = (194.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.05258,'amu*angstrom^2'), symmetry=1, barrier=(24.201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05258,'amu*angstrom^2'), symmetry=1, barrier=(24.201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33777,0.0398735,-3.70418e-05,1.81131e-08,-3.74289e-12,23480.8,16.0206], Tmin=(100,'K'), Tmax=(1118.07,'K')), NASAPolynomial(coeffs=[7.9573,0.0197688,-1.00692e-05,2.03004e-09,-1.46689e-13,22224.2,-11.7174], Tmin=(1118.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
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
    label = '[CH2]C([CH2])[CH][C]=O(2582)',
    structure = SMILES('[CH2]C([CH2])[CH][C]=O'),
    E0 = (477.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,318.813],'cm^-1')),
        HinderedRotor(inertia=(0.0281889,'amu*angstrom^2'), symmetry=1, barrier=(69.0787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00379562,'amu*angstrom^2'), symmetry=1, barrier=(9.32132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281264,'amu*angstrom^2'), symmetry=1, barrier=(69.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133159,'amu*angstrom^2'), symmetry=1, barrier=(9.3129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3452.89,'J/mol'), sigma=(6.01829,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.33 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41759,0.0538488,-5.61778e-05,3.36165e-08,-7.9844e-12,57553.1,26.3147], Tmin=(100,'K'), Tmax=(1127.64,'K')), NASAPolynomial(coeffs=[9.9417,0.0197888,-5.78525e-06,8.17611e-10,-4.62739e-14,55873.7,-14.7554], Tmin=(1127.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH][C]=O(9637)',
    structure = SMILES('[CH2][CH]C=[C][O]'),
    E0 = (502.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,458.276,458.348],'cm^-1')),
        HinderedRotor(inertia=(0.000816836,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000816843,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47043,0.0263284,1.39799e-06,-2.1241e-08,1.01461e-11,60468.5,21.3394], Tmin=(100,'K'), Tmax=(980.457,'K')), NASAPolynomial(coeffs=[9.92325,0.0131468,-4.7862e-06,8.81255e-10,-6.33441e-14,58179.2,-18.6908], Tmin=(980.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]([CH][C]=O)CO[O](10566)',
    structure = SMILES('[CH2][C]([CH][C]=O)CO[O]'),
    E0 = (527.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,2327.95,2328.01],'cm^-1')),
        HinderedRotor(inertia=(0.137564,'amu*angstrom^2'), symmetry=1, barrier=(3.16288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.94911,'amu*angstrom^2'), symmetry=1, barrier=(67.8059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0176322,'amu*angstrom^2'), symmetry=1, barrier=(67.807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137674,'amu*angstrom^2'), symmetry=1, barrier=(3.1654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137567,'amu*angstrom^2'), symmetry=1, barrier=(3.16295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01066,0.069968,-9.49528e-05,7.51799e-08,-2.38125e-11,63516,36.4543], Tmin=(100,'K'), Tmax=(879.576,'K')), NASAPolynomial(coeffs=[8.20702,0.0293151,-1.21072e-05,2.14238e-09,-1.41167e-13,62556.6,4.40245], Tmin=(879.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH][C]=O)[CH]O[O](10567)',
    structure = SMILES('[CH2]C([CH][C]=O)[CH]O[O]'),
    E0 = (524.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,1900.68,1900.73],'cm^-1')),
        HinderedRotor(inertia=(0.20578,'amu*angstrom^2'), symmetry=1, barrier=(4.73128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205726,'amu*angstrom^2'), symmetry=1, barrier=(4.73005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.83057,'amu*angstrom^2'), symmetry=1, barrier=(65.0804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205806,'amu*angstrom^2'), symmetry=1, barrier=(4.73188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.83065,'amu*angstrom^2'), symmetry=1, barrier=(65.0822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.408852,0.0862342,-0.000135361,1.1357e-07,-3.66044e-11,63254.8,33.2253], Tmin=(100,'K'), Tmax=(906.447,'K')), NASAPolynomial(coeffs=[9.57423,0.0287382,-1.2e-05,2.08938e-09,-1.34419e-13,62293.7,-6.22823], Tmin=(906.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](C[C]=O)CO[O](2912)',
    structure = SMILES('[CH2][C](C[C]=O)CO[O]'),
    E0 = (359.711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,180,2587.1],'cm^-1')),
        HinderedRotor(inertia=(0.0751877,'amu*angstrom^2'), symmetry=1, barrier=(1.72871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752082,'amu*angstrom^2'), symmetry=1, barrier=(1.72919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752012,'amu*angstrom^2'), symmetry=1, barrier=(1.72902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752177,'amu*angstrom^2'), symmetry=1, barrier=(1.7294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751839,'amu*angstrom^2'), symmetry=1, barrier=(1.72863,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869579,0.0760589,-0.000111481,9.64993e-08,-3.2996e-11,43369,36.5064], Tmin=(100,'K'), Tmax=(864.849,'K')), NASAPolynomial(coeffs=[6.36915,0.035763,-1.58179e-05,2.89067e-09,-1.93762e-13,42973.5,13.9856], Tmin=(864.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C[C]([CH][C]=O)CO[O](10568)',
    structure = SMILES('C[C]([CH][C]=O)CO[O]'),
    E0 = (322.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16579,0.0661189,-7.43407e-05,5.1472e-08,-1.51545e-11,38845.3,33.4314], Tmin=(100,'K'), Tmax=(813.358,'K')), NASAPolynomial(coeffs=[7.71266,0.0339209,-1.49588e-05,2.79796e-09,-1.93137e-13,37780.3,3.19929], Tmin=(813.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH][C]=O)[CH]OO(10569)',
    structure = SMILES('[CH2]C([CH][C]=O)[CH]OO'),
    E0 = (372.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16439,0.08958,-0.000129728,1.02996e-07,-3.22964e-11,44983.5,33.345], Tmin=(100,'K'), Tmax=(870.283,'K')), NASAPolynomial(coeffs=[11.2095,0.0303188,-1.29449e-05,2.3192e-09,-1.5364e-13,43382.8,-16.5586], Tmin=(870.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOOH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]O[O])C[C]=O(2914)',
    structure = SMILES('[CH2]C([CH]O[O])C[C]=O'),
    E0 = (357.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,773.561],'cm^-1')),
        HinderedRotor(inertia=(0.198551,'amu*angstrom^2'), symmetry=1, barrier=(4.56508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198726,'amu*angstrom^2'), symmetry=1, barrier=(4.56911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199025,'amu*angstrom^2'), symmetry=1, barrier=(4.57597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198902,'amu*angstrom^2'), symmetry=1, barrier=(4.57315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198459,'amu*angstrom^2'), symmetry=1, barrier=(4.56296,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269904,0.092293,-0.000151735,1.34611e-07,-4.56231e-11,43107.8,33.9632], Tmin=(100,'K'), Tmax=(888.454,'K')), NASAPolynomial(coeffs=[7.75455,0.0351552,-1.56929e-05,2.83347e-09,-1.86665e-13,42702.9,3.94578], Tmin=(888.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([CH][C]=O)[CH]O[O](10570)',
    structure = SMILES('CC([CH][C]=O)[CH]O[O]'),
    E0 = (319.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450807,0.0838431,-0.000120457,9.81997e-08,-3.19503e-11,38588.9,31.2937], Tmin=(100,'K'), Tmax=(848.354,'K')), NASAPolynomial(coeffs=[9.41501,0.032745,-1.44936e-05,2.65818e-09,-1.7905e-13,37385.7,-8.60607], Tmin=(848.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJCHO) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([CH]C=O)[CH]O[O](10572)',
    structure = SMILES('[CH2]C([CH]O[O])C=C[O]'),
    E0 = (341.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,293.047,293.054,3046.48],'cm^-1')),
        HinderedRotor(inertia=(0.287272,'amu*angstrom^2'), symmetry=1, barrier=(17.5081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27045,'amu*angstrom^2'), symmetry=1, barrier=(77.415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287249,'amu*angstrom^2'), symmetry=1, barrier=(17.5077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2704,'amu*angstrom^2'), symmetry=1, barrier=(77.4152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.135225,0.0802326,-8.90688e-05,5.11517e-08,-1.14954e-11,41245.5,33.3721], Tmin=(100,'K'), Tmax=(1094.63,'K')), NASAPolynomial(coeffs=[16.5331,0.0203108,-6.95531e-06,1.14123e-09,-7.34796e-14,37655.6,-47.2205], Tmin=(1094.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH][C]=O)COO(10573)',
    structure = SMILES('[CH2][C]([CH][C]=O)COO'),
    E0 = (375.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180.14,1603.9],'cm^-1')),
        HinderedRotor(inertia=(0.163169,'amu*angstrom^2'), symmetry=1, barrier=(3.75641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163169,'amu*angstrom^2'), symmetry=1, barrier=(3.75641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163169,'amu*angstrom^2'), symmetry=1, barrier=(3.75641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163169,'amu*angstrom^2'), symmetry=1, barrier=(3.75641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163169,'amu*angstrom^2'), symmetry=1, barrier=(3.75641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163169,'amu*angstrom^2'), symmetry=1, barrier=(3.75641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856945,0.072157,-8.4849e-05,5.81563e-08,-1.64424e-11,45240.8,36.2541], Tmin=(100,'K'), Tmax=(855.83,'K')), NASAPolynomial(coeffs=[9.58867,0.0313462,-1.33201e-05,2.43692e-09,-1.65847e-13,43746.3,-4.51202], Tmin=(855.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CJCOOH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(C=C=O)CO[O](2911)',
    structure = SMILES('[CH2]C(C=C=O)CO[O]'),
    E0 = (139.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.309189,'amu*angstrom^2'), symmetry=1, barrier=(7.10887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30918,'amu*angstrom^2'), symmetry=1, barrier=(7.10865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309165,'amu*angstrom^2'), symmetry=1, barrier=(7.10831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598436,'amu*angstrom^2'), symmetry=1, barrier=(13.7592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.188157,0.0918585,-0.000143748,1.23191e-07,-4.13894e-11,16954.6,30.9627], Tmin=(100,'K'), Tmax=(850.187,'K')), NASAPolynomial(coeffs=[9.5358,0.0335606,-1.56294e-05,2.91836e-09,-1.97821e-13,15882.6,-9.57379], Tmin=(850.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CJC(C)C=C=O)"""),
)

species(
    label = '[O]OC[CH]C[CH][C]=O(2999)',
    structure = SMILES('[O]OC[CH]C[CH][C]=O'),
    E0 = (340.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428438,0.0860431,-0.000131112,1.11937e-07,-3.75725e-11,41098.6,33.8085], Tmin=(100,'K'), Tmax=(857.884,'K')), NASAPolynomial(coeffs=[8.56584,0.0339855,-1.54087e-05,2.84325e-09,-1.91483e-13,40221.9,-1.17494], Tmin=(857.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCOOH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OCC1C=C([O])C1(10574)',
    structure = SMILES('[O]OCC1[CH]C(=O)C1'),
    E0 = (77.0597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28558,0.0482372,-5.97171e-06,-2.46093e-08,1.20294e-11,9375.81,27.0835], Tmin=(100,'K'), Tmax=(1030.88,'K')), NASAPolynomial(coeffs=[13.3504,0.0264643,-1.07267e-05,2.02887e-09,-1.45007e-13,5557.77,-37.9429], Tmin=(1030.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.0597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(ROOJ) + radical(CCJC=O)"""),
)

species(
    label = '[CH2]C1C=C([O])OOC1(10575)',
    structure = SMILES('[CH2]C1[CH]C(=O)OOC1'),
    E0 = (-57.8708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56153,0.0272829,8.22405e-05,-1.27778e-07,5.01135e-11,-6848.9,25.1154], Tmin=(100,'K'), Tmax=(983.091,'K')), NASAPolynomial(coeffs=[19.269,0.0209987,-8.51351e-06,1.81082e-09,-1.44978e-13,-13508.4,-76.1738], Tmin=(983.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.8708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C=C([O])OC1(10576)',
    structure = SMILES('[CH2]C1[CH]C(=O)OC1'),
    E0 = (-9.94173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53139,0.0144222,8.09207e-05,-1.18522e-07,4.85664e-11,-1126.11,24.2226], Tmin=(100,'K'), Tmax=(904.219,'K')), NASAPolynomial(coeffs=[12.1035,0.0185852,-3.13609e-06,3.33784e-10,-2.18389e-14,-4758.43,-31.507], Tmin=(904.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.94173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(butyrolactone) + radical(CCJCO) + radical(Isobutyl)"""),
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
    label = '[CH2]C(C#C[O])CO[O](10578)',
    structure = SMILES('[CH2]C([C]=C=O)CO[O]'),
    E0 = (345.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.348708,'amu*angstrom^2'), symmetry=1, barrier=(8.01749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34764,'amu*angstrom^2'), symmetry=1, barrier=(7.99292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348065,'amu*angstrom^2'), symmetry=1, barrier=(8.0027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.31758,'amu*angstrom^2'), symmetry=1, barrier=(76.2776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.131238,0.095918,-0.000164245,1.47314e-07,-5.04864e-11,41705,30.9969], Tmin=(100,'K'), Tmax=(869.811,'K')), NASAPolynomial(coeffs=[8.77816,0.0330014,-1.58179e-05,2.95148e-09,-1.98375e-13,41076.6,-4.47894], Tmin=(869.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CJC(C)C=C=O) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = '[CH2]C([C]=C[O])CO[O](10579)',
    structure = SMILES('[CH2]C([C]=C[O])CO[O]'),
    E0 = (391.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,380.731,384.048,384.976],'cm^-1')),
        HinderedRotor(inertia=(0.0974736,'amu*angstrom^2'), symmetry=1, barrier=(10.1967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.098774,'amu*angstrom^2'), symmetry=1, barrier=(10.2083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971754,'amu*angstrom^2'), symmetry=1, barrier=(10.2022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852109,'amu*angstrom^2'), symmetry=1, barrier=(88.3691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.170056,0.0812441,-9.40469e-05,5.65997e-08,-1.3346e-11,47167.2,33.5723], Tmin=(100,'K'), Tmax=(1043.47,'K')), NASAPolynomial(coeffs=[15.8986,0.0209512,-7.37552e-06,1.22606e-09,-7.93356e-14,43884.8,-42.9781], Tmin=(1043.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'CC([C]=[C][O])CO[O](10580)',
    structure = SMILES('CC([C][C]=O)CO[O]'),
    E0 = (411.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.246877,0.0903358,-0.000139989,1.19637e-07,-4.01833e-11,49671.8,32.0201], Tmin=(100,'K'), Tmax=(849.273,'K')), NASAPolynomial(coeffs=[9.3385,0.0337272,-1.56534e-05,2.91983e-09,-1.9791e-13,48624.7,-7.42894], Tmin=(849.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=[C]O)CO[O](10581)',
    structure = SMILES('[CH2]C([C]=[C]O)CO[O]'),
    E0 = (489.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,315.807,317.34],'cm^-1')),
        HinderedRotor(inertia=(0.114746,'amu*angstrom^2'), symmetry=1, barrier=(8.14062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114367,'amu*angstrom^2'), symmetry=1, barrier=(8.14821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114638,'amu*angstrom^2'), symmetry=1, barrier=(8.14732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115096,'amu*angstrom^2'), symmetry=1, barrier=(8.1458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32568,'amu*angstrom^2'), symmetry=1, barrier=(23.1986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190778,0.0869811,-0.000111844,6.93151e-08,-1.42805e-11,58981.7,35.2389], Tmin=(100,'K'), Tmax=(735.01,'K')), NASAPolynomial(coeffs=[14.7695,0.0218386,-7.87357e-06,1.29069e-09,-8.12547e-14,56455.1,-33.2154], Tmin=(735.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Isobutyl) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C=[C]O)CO[O](10582)',
    structure = SMILES('[CH2][C](C=[C]O)CO[O]'),
    E0 = (404.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,353.332,353.943],'cm^-1')),
        HinderedRotor(inertia=(0.0770051,'amu*angstrom^2'), symmetry=1, barrier=(6.82136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134636,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180559,'amu*angstrom^2'), symmetry=1, barrier=(15.9682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0767419,'amu*angstrom^2'), symmetry=1, barrier=(6.82108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767013,'amu*angstrom^2'), symmetry=1, barrier=(68.3588,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291513,0.0784892,-8.73029e-05,4.89031e-08,-1.01547e-11,48728.3,35.2708], Tmin=(100,'K'), Tmax=(913.542,'K')), NASAPolynomial(coeffs=[15.6701,0.0202983,-6.77194e-06,1.09295e-09,-6.97975e-14,45536.9,-39.6201], Tmin=(913.542,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([C]=[C][O])COO(10583)',
    structure = SMILES('[CH2]C([C][C]=O)COO'),
    E0 = (465.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0334483,0.095994,-0.000148949,1.23972e-07,-4.03024e-11,56066.1,34.05], Tmin=(100,'K'), Tmax=(865.275,'K')), NASAPolynomial(coeffs=[11.1212,0.0313224,-1.41175e-05,2.58395e-09,-1.72763e-13,54626.4,-15.3155], Tmin=(865.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]O[O])C=[C]O(10584)',
    structure = SMILES('[CH2]C([CH]O[O])C=[C]O'),
    E0 = (440.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,269.51,269.675],'cm^-1')),
        HinderedRotor(inertia=(0.242365,'amu*angstrom^2'), symmetry=1, barrier=(12.4578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243336,'amu*angstrom^2'), symmetry=1, barrier=(12.4635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242923,'amu*angstrom^2'), symmetry=1, barrier=(12.4594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49509,'amu*angstrom^2'), symmetry=1, barrier=(77.0298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.01034,'amu*angstrom^2'), symmetry=1, barrier=(77.0241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146504,0.0862426,-0.000108779,6.82306e-08,-1.54349e-11,53060.3,35.0636], Tmin=(100,'K'), Tmax=(787.498,'K')), NASAPolynomial(coeffs=[15.1347,0.0216459,-7.70709e-06,1.265e-09,-8.02522e-14,50342.1,-35.9357], Tmin=(787.498,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(C=CJO)"""),
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
    E0 = (336.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (534.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (582.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (581.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (837.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (774.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (828.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (791.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (341.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (343.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (343.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (419.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (419.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (475.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (359.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (359.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (361.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (445.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (500.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (497.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (433.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (469.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (708.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (739.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (736.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (517.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (477.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (502.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (500.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (453.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (489.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (438.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (459.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (336.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (586.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (344.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (344.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (391.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (399.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (572.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (544.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (595.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (563.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (639.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (487.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (506.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (466.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['HCCO(2227)', 'allylperoxy(462)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([C]=O)CO[O](10551)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OCC[CH][CH][C]=O(10552)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH2]C([CH][C]=O)C[O](10553)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[O]OC[CH][CH][C]=O(10554)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2]C([C][C]=O)CO[O](10555)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C([CH][C]=O)CO[O](10556)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[O]OCC1CC1[C]=O(10557)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[CH2]C1COOC1[C]=O(10558)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['O=[C][CH]C1COOC1(10559)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['O(4)', '[CH2]C1COC1[C]=O(10560)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['O(4)', 'O=[C][CH]C1COC1(10561)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['HO2(9)', '[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['C=C(C[C]=O)CO[O](2910)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['CC(=C[C]=O)CO[O](10562)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[CH2]C(=C[C]=O)COO(10563)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C(=C[C]=O)CO[O](10564)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(T)(20)', '[O]OCC=C[C]=O(10565)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C][O](6861)', 'allylperoxy(462)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(84.9359,'m^3/(mol*s)'), n=1.6115, Ea=(0.85772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;Y_1centerbirad] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]O[O](46)', 'C=CC=[C][O](9589)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(182.434,'m^3/(mol*s)'), n=0.88, Ea=(33.1163,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-OsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O2(2)', '[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.4339e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]O[O](46)', '[CH2][CH][CH][C]=O(9637)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][C]([CH][C]=O)CO[O](10566)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]C([CH][C]=O)[CH]O[O](10567)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C](C[C]=O)CO[O](2912)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['C[C]([CH][C]=O)CO[O](10568)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[CH2]C([CH][C]=O)[CH]OO(10569)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.66e+08,'s^-1'), n=1.28, Ea=(166.272,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 241 used for R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([CH]O[O])C[C]=O(2914)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['CC([CH][C]=O)[CH]O[O](10570)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.3e-15,'s^-1'), n=8.11, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 339 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[CH2][C]([CH]C=O)CO[O](10571)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([CH]C=O)[CH]O[O](10572)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeO;XH_out] for rate rule [R4HJ_2;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C]([CH][C]=O)COO(10573)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[CH2]C(C=C=O)CO[O](2911)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]OC[CH]C[CH][C]=O(2999)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[O]OCC1C=C([O])C1(10574)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['[CH2]C1C=C([O])OOC1(10575)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['O(4)', '[CH2]C1C=C([O])OC1(10576)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_DSS;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    products = ['C=C(C=C[O])CO[O](10577)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', '[CH2]C(C#C[O])CO[O](10578)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['HCCO(2227)', '[CH2][CH]CO[O](461)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([C]=C[O])CO[O](10579)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CC([C]=[C][O])CO[O](10580)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([C]=[C]O)CO[O](10581)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C](C=[C]O)CO[O](10582)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([C]=[C][O])COO(10583)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C([CH]O[O])C=[C]O(10584)'],
    products = ['[CH2]C([CH][C]=O)CO[O](2913)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_1H;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2349',
    isomers = [
        '[CH2]C([CH][C]=O)CO[O](2913)',
    ],
    reactants = [
        ('HCCO(2227)', 'allylperoxy(462)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2349',
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


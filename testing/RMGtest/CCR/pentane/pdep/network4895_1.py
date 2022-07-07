species(
    label = '[CH]=[C]OC([C]=C)O[O](21186)',
    structure = SMILES('[CH]=[C]OC([C]=C)O[O]'),
    E0 = (665.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,242.768,242.769,242.771],'cm^-1')),
        HinderedRotor(inertia=(0.365522,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365527,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365533,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365534,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330871,0.0844002,-0.000119246,8.53778e-08,-2.39875e-11,80148.5,32.599], Tmin=(100,'K'), Tmax=(875.533,'K')), NASAPolynomial(coeffs=[14.3528,0.0203398,-9.49635e-06,1.81049e-09,-1.25905e-13,77693.1,-33.1847], Tmin=(875.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OC=C=C(16179)',
    structure = SMILES('C#CO[CH][C]=C'),
    E0 = (482.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100,279.596,279.598,279.634],'cm^-1')),
        HinderedRotor(inertia=(0.594372,'amu*angstrom^2'), symmetry=1, barrier=(32.9751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594396,'amu*angstrom^2'), symmetry=1, barrier=(32.9753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594461,'amu*angstrom^2'), symmetry=1, barrier=(32.9754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36759,0.0525792,-5.21763e-05,2.56965e-08,-4.95256e-12,58089.6,20.1027], Tmin=(100,'K'), Tmax=(1266.22,'K')), NASAPolynomial(coeffs=[13.9321,0.0128878,-5.15679e-06,9.40687e-10,-6.48368e-14,54907.7,-43.4797], Tmin=(1266.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
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
    label = 'C=[C]C([O])O[O](21017)',
    structure = SMILES('C=[C]C([O])O[O]'),
    E0 = (333.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1685,370,180,2080.84],'cm^-1')),
        HinderedRotor(inertia=(0.384905,'amu*angstrom^2'), symmetry=1, barrier=(8.84973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383233,'amu*angstrom^2'), symmetry=1, barrier=(8.81127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88098,0.0575616,-0.000110039,1.12298e-07,-4.25317e-11,40131.1,23.4867], Tmin=(100,'K'), Tmax=(855.741,'K')), NASAPolynomial(coeffs=[2.29877,0.0295719,-1.53376e-05,2.96524e-09,-2.03584e-13,41012.9,27.1064], Tmin=(855.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_S)"""),
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
    label = '[CH]=[C]OC([O])[C]=C(23352)',
    structure = SMILES('[CH]=[C]OC([O])[C]=C'),
    E0 = (667.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,269.133,269.143,591.526,1137.75],'cm^-1')),
        HinderedRotor(inertia=(0.279916,'amu*angstrom^2'), symmetry=1, barrier=(14.3815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27973,'amu*angstrom^2'), symmetry=1, barrier=(14.3822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279833,'amu*angstrom^2'), symmetry=1, barrier=(14.382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11321,0.0682038,-9.5221e-05,7.21134e-08,-2.20429e-11,80383.3,29.1317], Tmin=(100,'K'), Tmax=(797.975,'K')), NASAPolynomial(coeffs=[9.98609,0.0237265,-1.1614e-05,2.26364e-09,-1.59321e-13,78967.3,-11.6724], Tmin=(797.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]O[CH]O[O](23353)',
    structure = SMILES('[CH]=[C]O[CH]O[O]'),
    E0 = (559.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.785168,'amu*angstrom^2'), symmetry=1, barrier=(18.0526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784505,'amu*angstrom^2'), symmetry=1, barrier=(18.0373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785235,'amu*angstrom^2'), symmetry=1, barrier=(18.0541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.968306,0.0647455,-9.59315e-05,6.52197e-08,-1.66106e-11,67390.5,24.4858], Tmin=(100,'K'), Tmax=(1035.87,'K')), NASAPolynomial(coeffs=[16.1825,0.00220607,1.17552e-07,-1.27832e-10,1.30247e-14,64441.8,-48.4686], Tmin=(1035.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(OCJO) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[C]=[C]OC([C]=C)O[O](23354)',
    structure = SMILES('[C]=[C]OC([C]=C)O[O]'),
    E0 = (976.327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1670,1700,300,440,314.016,314.022,314.234,314.918],'cm^-1')),
        HinderedRotor(inertia=(0.198994,'amu*angstrom^2'), symmetry=1, barrier=(13.941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199795,'amu*angstrom^2'), symmetry=1, barrier=(13.925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199721,'amu*angstrom^2'), symmetry=1, barrier=(13.9186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200919,'amu*angstrom^2'), symmetry=1, barrier=(13.925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.237383,0.0889256,-0.000146295,1.20954e-07,-3.88685e-11,117555,32.9009], Tmin=(100,'K'), Tmax=(832.521,'K')), NASAPolynomial(coeffs=[13.0898,0.0202605,-1.01208e-05,1.93385e-09,-1.32405e-13,115654,-25.31], Tmin=(832.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(976.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C1OC(O[O])C1=C(23355)',
    structure = SMILES('[CH]=C1OC(O[O])C1=C'),
    E0 = (306.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.54341,0.0673396,-6.23745e-05,2.51357e-08,-2.97996e-12,37014.4,23.0499], Tmin=(100,'K'), Tmax=(1053.15,'K')), NASAPolynomial(coeffs=[17.958,0.0143393,-5.60491e-06,1.04865e-09,-7.48889e-14,32617.6,-65.3276], Tmin=(1053.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC1OOC1=C(23356)',
    structure = SMILES('[CH]=[C]OC1OOC1=C'),
    E0 = (441.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.90686,0.052984,-1.41541e-05,-2.76549e-08,1.58812e-11,53266.7,28.3975], Tmin=(100,'K'), Tmax=(994.605,'K')), NASAPolynomial(coeffs=[19.5243,0.0125362,-5.07244e-06,1.05836e-09,-8.3346e-14,47860.5,-69.8804], Tmin=(994.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OOC([C]=C)O1(23357)',
    structure = SMILES('[CH]=C1OOC([C]=C)O1'),
    E0 = (415.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10062,0.0518539,-2.11467e-05,-1.35083e-08,9.36742e-12,50099.9,25.6631], Tmin=(100,'K'), Tmax=(1027.21,'K')), NASAPolynomial(coeffs=[16.5479,0.0170087,-7.21841e-06,1.43596e-09,-1.06742e-13,45591.2,-55.7747], Tmin=(1027.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC1OC1=C(23358)',
    structure = SMILES('[CH]=[C]OC1OC1=C'),
    E0 = (375.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801594,0.0520979,-7.91125e-06,-4.81492e-08,2.8059e-11,45266.3,22.072], Tmin=(100,'K'), Tmax=(922.121,'K')), NASAPolynomial(coeffs=[24.4337,-0.00282407,4.01533e-06,-8.03634e-10,4.92165e-14,38884.6,-100.994], Tmin=(922.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OC([C]=C)O1(23359)',
    structure = SMILES('[CH]=C1OC([C]=C)O1'),
    E0 = (391.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1273,0.0457026,2.89563e-06,-5.00539e-08,2.57506e-11,47221.8,20.7091], Tmin=(100,'K'), Tmax=(951.234,'K')), NASAPolynomial(coeffs=[21.2377,0.00445757,-3.77202e-07,1.1597e-10,-1.75261e-14,41435.9,-85.6086], Tmin=(951.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=COC(=C=C)O[O](23360)',
    structure = SMILES('[CH]=COC(=C=C)O[O]'),
    E0 = (476.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444299,0.0855476,-0.000139031,1.18916e-07,-3.95456e-11,57473,29.6535], Tmin=(100,'K'), Tmax=(851.667,'K')), NASAPolynomial(coeffs=[10.2704,0.0266057,-1.26899e-05,2.38338e-09,-1.61654e-13,56263.2,-13.4507], Tmin=(851.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(=C=C)OO(23361)',
    structure = SMILES('C#CO[C]([C]=C)OO'),
    E0 = (441.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0422973,0.0882158,-0.00011453,7.18838e-08,-1.75323e-11,53187.6,27.4882], Tmin=(100,'K'), Tmax=(1008.21,'K')), NASAPolynomial(coeffs=[18.036,0.0168291,-8.32482e-06,1.65888e-09,-1.19525e-13,49559.2,-59.4685], Tmin=(1008.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cs_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]OC(=C=C)O[O](23362)',
    structure = SMILES('C#CO[C]([C]=C)O[O]'),
    E0 = (593.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.964518,'amu*angstrom^2'), symmetry=1, barrier=(22.1762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964724,'amu*angstrom^2'), symmetry=1, barrier=(22.1809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964883,'amu*angstrom^2'), symmetry=1, barrier=(22.1846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.964646,'amu*angstrom^2'), symmetry=1, barrier=(22.1791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.250798,0.0852754,-0.000121522,8.42201e-08,-2.26311e-11,71460.5,27.499], Tmin=(100,'K'), Tmax=(918.376,'K')), NASAPolynomial(coeffs=[16.2767,0.0154751,-7.51755e-06,1.46286e-09,-1.0319e-13,68516.9,-48.4524], Tmin=(918.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(Cs_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]OC(C#C)O[O](23363)',
    structure = SMILES('[CH]=[C]OC(C#C)O[O]'),
    E0 = (593.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2175,525,1685,370,750,770,3400,2100,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02928,'amu*angstrom^2'), symmetry=1, barrier=(23.6652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02911,'amu*angstrom^2'), symmetry=1, barrier=(23.6612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02921,'amu*angstrom^2'), symmetry=1, barrier=(23.6635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02929,'amu*angstrom^2'), symmetry=1, barrier=(23.6654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194887,0.0844408,-0.00011984,8.15208e-08,-2.08617e-11,71509.5,30.2254], Tmin=(100,'K'), Tmax=(838.56,'K')), NASAPolynomial(coeffs=[16.9644,0.012172,-4.38145e-06,7.12668e-10,-4.46328e-14,68425.5,-49.344], Tmin=(838.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]O[CH][C]=C(16181)',
    structure = SMILES('[CH]=[C]O[CH][C]=C'),
    E0 = (759.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,180,597.029,833.203],'cm^-1')),
        HinderedRotor(inertia=(0.74415,'amu*angstrom^2'), symmetry=1, barrier=(17.1095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743464,'amu*angstrom^2'), symmetry=1, barrier=(17.0937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191273,'amu*angstrom^2'), symmetry=1, barrier=(95.2621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3208.06,'J/mol'), sigma=(5.53191,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.09 K, Pc=43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69949,0.0499227,-5.16581e-05,2.80139e-08,-6.09932e-12,91451,24.9866], Tmin=(100,'K'), Tmax=(1109.68,'K')), NASAPolynomial(coeffs=[10.8034,0.0171054,-7.2963e-06,1.3616e-09,-9.46263e-14,89430.6,-19.8815], Tmin=(1109.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(759.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C]([C]=C)O[O](23364)',
    structure = SMILES('[CH]=[C]O[C]([C]=C)O[O]'),
    E0 = (870.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,322.234,322.268,322.382],'cm^-1')),
        HinderedRotor(inertia=(0.192533,'amu*angstrom^2'), symmetry=1, barrier=(14.2013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192808,'amu*angstrom^2'), symmetry=1, barrier=(14.2013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192546,'amu*angstrom^2'), symmetry=1, barrier=(14.2027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192811,'amu*angstrom^2'), symmetry=1, barrier=(14.2017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.287447,0.086151,-0.000133583,1.03182e-07,-3.1043e-11,104835,33.4384], Tmin=(100,'K'), Tmax=(820.746,'K')), NASAPolynomial(coeffs=[14.33,0.0177212,-8.53583e-06,1.62161e-09,-1.11402e-13,102529,-31.5366], Tmin=(820.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(870.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cs_P) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC([C]=[CH])O[O](23365)',
    structure = SMILES('[CH]=[C]OC([C]=[CH])O[O]'),
    E0 = (912.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,1670,1700,300,440,3115,3125,620,680,785,800,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.693125,'amu*angstrom^2'), symmetry=1, barrier=(15.9363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693002,'amu*angstrom^2'), symmetry=1, barrier=(15.9335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.692981,'amu*angstrom^2'), symmetry=1, barrier=(15.933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693079,'amu*angstrom^2'), symmetry=1, barrier=(15.9352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.147353,0.0893071,-0.000141187,1.10291e-07,-3.3444e-11,109873,33.6548], Tmin=(100,'K'), Tmax=(815.444,'K')), NASAPolynomial(coeffs=[14.8937,0.0169781,-8.15e-06,1.53666e-09,-1.0473e-13,107468,-34.4805], Tmin=(815.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(912.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C](C=C)O[O](23366)',
    structure = SMILES('[CH]=[C]O[C](C=C)O[O]'),
    E0 = (632.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.417471,0.0803606,-0.00010438,6.77358e-08,-1.71967e-11,76227.1,32.501], Tmin=(100,'K'), Tmax=(967.82,'K')), NASAPolynomial(coeffs=[15.4053,0.0184144,-8.36894e-06,1.5986e-09,-1.12218e-13,73326,-39.316], Tmin=(967.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cs_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(C=[CH])O[O](23367)',
    structure = SMILES('[CH]=[C]OC(C=[CH])O[O]'),
    E0 = (674.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.732916,'amu*angstrom^2'), symmetry=1, barrier=(16.8512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732888,'amu*angstrom^2'), symmetry=1, barrier=(16.8505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732912,'amu*angstrom^2'), symmetry=1, barrier=(16.8511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732899,'amu*angstrom^2'), symmetry=1, barrier=(16.8508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289524,0.0833497,-0.000111285,7.37758e-08,-1.9067e-11,81264.9,32.6755], Tmin=(100,'K'), Tmax=(952.616,'K')), NASAPolynomial(coeffs=[15.9017,0.0177929,-8.05634e-06,1.53151e-09,-1.07063e-13,78290.5,-41.8862], Tmin=(952.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CO[C]([C]=C)O[O](23368)',
    structure = SMILES('[CH]=CO[C]([C]=C)O[O]'),
    E0 = (630.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.040279,0.0848686,-0.000108483,6.67603e-08,-1.57957e-11,76015.1,31.1941], Tmin=(100,'K'), Tmax=(1045.66,'K')), NASAPolynomial(coeffs=[18.8425,0.0129429,-5.30411e-06,9.77447e-10,-6.79111e-14,72083,-60.355], Tmin=(1045.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(630.823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cs_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C]([C]=C)OO(23369)',
    structure = SMILES('[CH]=[C]O[C]([C]=C)OO'),
    E0 = (718.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.175059,0.0878436,-0.000121688,8.37101e-08,-2.25488e-11,86557.7,33.0904], Tmin=(100,'K'), Tmax=(912.335,'K')), NASAPolynomial(coeffs=[15.705,0.0197571,-9.74848e-06,1.91554e-09,-1.35995e-13,83723.9,-40.4082], Tmin=(912.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(718.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cs_P) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]O[C]([C]=C)O[O](21810)',
    structure = SMILES('C=[C]O[C]([C]=C)O[O]'),
    E0 = (623.471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,260.816,260.817,260.818,260.819],'cm^-1')),
        HinderedRotor(inertia=(0.00247814,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333418,'amu*angstrom^2'), symmetry=1, barrier=(16.0947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333415,'amu*angstrom^2'), symmetry=1, barrier=(16.0947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333403,'amu*angstrom^2'), symmetry=1, barrier=(16.0947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467554,0.0812961,-0.000111884,7.8669e-08,-2.17995e-11,75110.3,32.3941], Tmin=(100,'K'), Tmax=(885.857,'K')), NASAPolynomial(coeffs=[13.8131,0.0210383,-9.85477e-06,1.88867e-09,-1.31996e-13,72745.8,-30.3734], Tmin=(885.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cs_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]C(O[O])OC=[CH](23370)',
    structure = SMILES('[CH]=[C]C(O[O])OC=[CH]'),
    E0 = (672.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1685,370,3115,3125,620,680,785,800,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.891185,'amu*angstrom^2'), symmetry=1, barrier=(20.4901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889582,'amu*angstrom^2'), symmetry=1, barrier=(20.4532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888255,'amu*angstrom^2'), symmetry=1, barrier=(20.4227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891039,'amu*angstrom^2'), symmetry=1, barrier=(20.4867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0746169,0.0877017,-0.000114836,7.20788e-08,-1.73571e-11,81052.4,31.3219], Tmin=(100,'K'), Tmax=(1028.68,'K')), NASAPolynomial(coeffs=[19.2894,0.0124049,-5.03936e-06,9.21599e-10,-6.36853e-14,77068.5,-62.6457], Tmin=(1028.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC([C]=[CH])OO(23371)',
    structure = SMILES('[CH]=[C]OC([C]=[CH])OO'),
    E0 = (760.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,1670,1700,300,440,3115,3125,620,680,785,800,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.899044,'amu*angstrom^2'), symmetry=1, barrier=(20.6708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899045,'amu*angstrom^2'), symmetry=1, barrier=(20.6708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899035,'amu*angstrom^2'), symmetry=1, barrier=(20.6706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899021,'amu*angstrom^2'), symmetry=1, barrier=(20.6703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.899011,'amu*angstrom^2'), symmetry=1, barrier=(20.67,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0404434,0.0909184,-0.000128925,9.02223e-08,-2.46385e-11,91595.8,33.2883], Tmin=(100,'K'), Tmax=(901.587,'K')), NASAPolynomial(coeffs=[16.2311,0.0190833,-9.40497e-06,1.84101e-09,-1.30215e-13,88676.5,-43.1442], Tmin=(901.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(O[O])O[C]=C(21811)',
    structure = SMILES('[CH]=[C]C(O[O])O[C]=C'),
    E0 = (665.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,242.768,242.769,242.771],'cm^-1')),
        HinderedRotor(inertia=(0.365522,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365527,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365533,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365534,'amu*angstrom^2'), symmetry=1, barrier=(15.2882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330871,0.0844002,-0.000119246,8.53778e-08,-2.39875e-11,80148.5,32.599], Tmin=(100,'K'), Tmax=(875.533,'K')), NASAPolynomial(coeffs=[14.3528,0.0203398,-9.49635e-06,1.81049e-09,-1.25905e-13,77693.1,-33.1847], Tmin=(875.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (665.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (665.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1180.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (910.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1215.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1188.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (673.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (673.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (672.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (796.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (772.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (743.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (743.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (817.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (818.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (665.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (665.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (862.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (751.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1082.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1124.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (847.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (780.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (818.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (853.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (874.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (705.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (793.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (856.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['HCCO(2227)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['O2(S)(666)', '[CH]=[C]OC=C=C(16179)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=[CH](9646)', 'C=[C]C([O])O[O](21017)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]OC([O])[C]=C(23352)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H2CC(T)(1341)', '[CH]=[C]O[CH]O[O](23353)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/O2;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[C]=[C]OC([C]=C)O[O](23354)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=C1OC(O[O])C1=C(23355)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=[C]OC1OOC1=C(23356)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=C1OOC([C]=C)O1(23357)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['O(4)', '[CH]=[C]OC1OC1=C(23358)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO_S;Y_rad_intra;OOJ] for rate rule [R2OO_S;Cd_rad_out;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['O(4)', '[CH]=C1OC([C]=C)O1(23359)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.67697e+10,'s^-1'), n=0.54, Ea=(107.243,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R3OO;Y_rad_intra;OO] for rate rule [R3OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=COC(=C=C)O[O](23360)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=[C]OC(=C=C)OO(23361)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]=[C]OC(=C=C)O[O](23362)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH]=[C]OC(C#C)O[O](23363)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HCCO(2227)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(62.8203,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 58.3 to 62.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[CH]=[C]OC=C=C(16179)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.973,'m^3/(mol*s)'), n=2.037, Ea=(191.783,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;OJ] for rate rule [Cds_Ca;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 191.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C][O](6861)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O2(2)', '[CH]=[C]O[CH][C]=C(16181)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH]=[C]O[C]([C]=C)O[O](23364)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH]=[C]OC([C]=[CH])O[O](23365)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=[C]O[C](C=C)O[O](23366)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]OC(C=[CH])O[O](23367)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['[CH]=CO[C]([C]=C)O[O](23368)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.38739e+08,'s^-1'), n=1.16185, Ea=(153.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]O[C]([C]=C)OO(23369)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    products = ['C=[C]O[C]([C]=C)O[O](21810)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]C(O[O])OC=[CH](23370)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]OC([C]=[CH])OO(23371)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]C(O[O])O[C]=C(21811)'],
    products = ['[CH]=[C]OC([C]=C)O[O](21186)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6Hall;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4895',
    isomers = [
        '[CH]=[C]OC([C]=C)O[O](21186)',
    ],
    reactants = [
        ('HCCO(2227)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4895',
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


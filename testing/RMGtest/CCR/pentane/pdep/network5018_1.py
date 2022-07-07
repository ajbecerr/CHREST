species(
    label = '[CH2][C]=CC([C]=C)O[O](21296)',
    structure = SMILES('C=[C][CH]C([C]=C)O[O]'),
    E0 = (660.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,376.897,376.905,376.921],'cm^-1')),
        HinderedRotor(inertia=(0.00118675,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27058,'amu*angstrom^2'), symmetry=1, barrier=(27.2787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270611,'amu*angstrom^2'), symmetry=1, barrier=(27.2793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270563,'amu*angstrom^2'), symmetry=1, barrier=(27.2789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366681,0.06927,-5.96723e-05,2.54595e-08,-4.27754e-12,79547,34.6298], Tmin=(100,'K'), Tmax=(1437.8,'K')), NASAPolynomial(coeffs=[18.2373,0.0195533,-7.80471e-06,1.41e-09,-9.58787e-14,74408.1,-58.0749], Tmin=(1437.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S)"""),
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
    label = '[CH2][C]=CC=C=C(16910)',
    structure = SMILES('C=[C]C=C[C]=C'),
    E0 = (543.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62996,'amu*angstrom^2'), symmetry=1, barrier=(37.4761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62989,'amu*angstrom^2'), symmetry=1, barrier=(37.4743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759224,0.0600512,-5.80613e-05,2.91517e-08,-5.61572e-12,65516.4,20.7626], Tmin=(100,'K'), Tmax=(1428.68,'K')), NASAPolynomial(coeffs=[14.7312,0.0143871,-3.24523e-06,3.6595e-10,-1.74458e-14,62192.1,-49.2908], Tmin=(1428.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = 'C=[C][CH]C([O])[C]=C(22553)',
    structure = SMILES('C=[C][CH]C([O])[C]=C'),
    E0 = (667.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,563.253,563.333,563.454,563.546],'cm^-1')),
        HinderedRotor(inertia=(0.112792,'amu*angstrom^2'), symmetry=1, barrier=(25.4289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112652,'amu*angstrom^2'), symmetry=1, barrier=(25.4196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112561,'amu*angstrom^2'), symmetry=1, barrier=(25.4234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17164,0.0494395,-1.02258e-05,-2.5621e-08,1.39273e-11,80344.5,30.0144], Tmin=(100,'K'), Tmax=(998.569,'K')), NASAPolynomial(coeffs=[16.1255,0.0184835,-7.20527e-06,1.39072e-09,-1.02721e-13,75914.8,-49.3339], Tmin=(998.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S)"""),
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
    label = 'C=[C][CH][CH]O[O](22554)',
    structure = SMILES('[CH2][C]=C[CH]O[O]'),
    E0 = (547.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,376.8],'cm^-1')),
        HinderedRotor(inertia=(0.00118735,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343266,'amu*angstrom^2'), symmetry=1, barrier=(34.5828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529098,'amu*angstrom^2'), symmetry=1, barrier=(53.3043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62471,0.0456463,-3.90107e-05,1.68223e-08,-2.87474e-12,65892.9,24.1318], Tmin=(100,'K'), Tmax=(1410.65,'K')), NASAPolynomial(coeffs=[12.7508,0.0140978,-5.46438e-06,9.68594e-10,-6.51325e-14,62753.9,-33.3731], Tmin=(1410.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = 'C=[C][C]C([C]=C)O[O](22555)',
    structure = SMILES('[CH2][C]=[C]C([C]=C)O[O]'),
    E0 = (963.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,303.573,303.66],'cm^-1')),
        HinderedRotor(inertia=(0.126206,'amu*angstrom^2'), symmetry=1, barrier=(8.25408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126133,'amu*angstrom^2'), symmetry=1, barrier=(8.25387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126143,'amu*angstrom^2'), symmetry=1, barrier=(8.2547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.450634,'amu*angstrom^2'), symmetry=1, barrier=(29.4811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503145,0.0838429,-0.000128076,1.09193e-07,-3.71979e-11,115993,33.4807], Tmin=(100,'K'), Tmax=(809.027,'K')), NASAPolynomial(coeffs=[9.1988,0.0320103,-1.55855e-05,2.99194e-09,-2.07105e-13,114875,-4.84022], Tmin=(809.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(963.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C(=C)C1O[O](22556)',
    structure = SMILES('C=[C]C1C(=C)C1O[O]'),
    E0 = (547.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787709,0.0693833,-6.63414e-05,3.34422e-08,-6.83525e-12,65928.6,25.2222], Tmin=(100,'K'), Tmax=(1170.65,'K')), NASAPolynomial(coeffs=[13.3463,0.0264721,-1.1358e-05,2.13012e-09,-1.48385e-13,62988.2,-37.3443], Tmin=(1170.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1OOC1[C]=C(22557)',
    structure = SMILES('C=[C]C1OOC1[C]=C'),
    E0 = (594.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928562,0.0613428,-4.66003e-05,1.75028e-08,-2.64118e-12,71612.2,26.7295], Tmin=(100,'K'), Tmax=(1552.33,'K')), NASAPolynomial(coeffs=[15.499,0.0237978,-1.03206e-05,1.92191e-09,-1.31889e-13,67088.6,-49.9718], Tmin=(1552.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(594.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C1OOC1=C(22558)',
    structure = SMILES('C=[C][CH]C1OOC1=C'),
    E0 = (436.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41225,0.0324167,6.38718e-05,-1.1051e-07,4.48817e-11,52644.8,28.738], Tmin=(100,'K'), Tmax=(976.184,'K')), NASAPolynomial(coeffs=[19.8872,0.0174965,-6.6002e-06,1.40258e-09,-1.14132e-13,46141.7,-74.7815], Tmin=(976.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1OC1[C]=C(22559)',
    structure = SMILES('C=[C]C1OC1[C]=C'),
    E0 = (543.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40093,0.0458208,-3.34401e-07,-4.26506e-08,2.38639e-11,65431.6,24.8067], Tmin=(100,'K'), Tmax=(882.601,'K')), NASAPolynomial(coeffs=[15.4565,0.0144536,-1.97624e-06,9.66594e-11,-1.45121e-15,61691.1,-48.3822], Tmin=(882.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C1OC1=C(22560)',
    structure = SMILES('C=[C][CH]C1OC1=C'),
    E0 = (370.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28851,0.0317549,6.93009e-05,-1.29907e-07,5.65726e-11,44645.2,22.4783], Tmin=(100,'K'), Tmax=(933.674,'K')), NASAPolynomial(coeffs=[24.8721,0.00200816,2.56134e-06,-4.76826e-10,1.98742e-14,37134,-106.32], Tmin=(933.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(=C=C)O[O](22300)',
    structure = SMILES('C=[C]CC(=C=C)O[O]'),
    E0 = (529.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.369063,'amu*angstrom^2'), symmetry=1, barrier=(8.48549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368691,'amu*angstrom^2'), symmetry=1, barrier=(8.47693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367804,'amu*angstrom^2'), symmetry=1, barrier=(8.45653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788274,0.0756188,-9.6937e-05,7.24186e-08,-2.23784e-11,63748.4,30.5538], Tmin=(100,'K'), Tmax=(783.664,'K')), NASAPolynomial(coeffs=[9.15179,0.0329291,-1.52245e-05,2.90478e-09,-2.02351e-13,62437.6,-7.75643], Tmin=(783.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C(C=C)O[O](22561)',
    structure = SMILES('[CH2]C=C(C=C=C)O[O]'),
    E0 = (379.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788401,0.0714444,-7.48016e-05,4.32747e-08,-1.02254e-11,45714.6,27.7965], Tmin=(100,'K'), Tmax=(1019.18,'K')), NASAPolynomial(coeffs=[11.6574,0.0287864,-1.20183e-05,2.20673e-09,-1.51533e-13,43499.1,-24.8466], Tmin=(1019.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C=C([C]=C)OO(22562)',
    structure = SMILES('C=[C]C=C([C]=C)OO'),
    E0 = (454.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0269535,0.0907666,-0.0001184,7.84059e-08,-1.91951e-11,54778.1,27.3861], Tmin=(100,'K'), Tmax=(755.135,'K')), NASAPolynomial(coeffs=[14.4262,0.025345,-1.00033e-05,1.74e-09,-1.14258e-13,52294,-40.0866], Tmin=(755.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=C([C]=C)O[O](22563)',
    structure = SMILES('C=[C]C=C([C]=C)O[O]'),
    E0 = (606.291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07567,'amu*angstrom^2'), symmetry=1, barrier=(24.7317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07573,'amu*angstrom^2'), symmetry=1, barrier=(24.7332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07568,'amu*angstrom^2'), symmetry=1, barrier=(24.7319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.11161,0.0895473,-0.000132775,1.02608e-07,-3.05779e-11,73056.2,27.8252], Tmin=(100,'K'), Tmax=(931.466,'K')), NASAPolynomial(coeffs=[13.0126,0.0233573,-8.81074e-06,1.4494e-09,-8.985e-14,71120.8,-30.9862], Tmin=(931.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC([CH][C]=C)O[O](22564)',
    structure = SMILES('C#CC([CH][C]=C)O[O]'),
    E0 = (635.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.47346,'amu*angstrom^2'), symmetry=1, barrier=(33.8777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47405,'amu*angstrom^2'), symmetry=1, barrier=(33.8913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47228,'amu*angstrom^2'), symmetry=1, barrier=(33.8507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4747,'amu*angstrom^2'), symmetry=1, barrier=(33.9064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286506,0.0813088,-9.71669e-05,5.98741e-08,-1.45378e-11,76520.4,28.4513], Tmin=(100,'K'), Tmax=(1009.36,'K')), NASAPolynomial(coeffs=[15.2621,0.0219626,-8.97401e-06,1.62476e-09,-1.10706e-13,73497.2,-43.9371], Tmin=(1009.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(635.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C([C]=C)O[O](22565)',
    structure = SMILES('[CH]=C=CC([C]=C)O[O]'),
    E0 = (667.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,540,610,2055,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.752132,'amu*angstrom^2'), symmetry=1, barrier=(17.293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.751626,'amu*angstrom^2'), symmetry=1, barrier=(17.2814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.751904,'amu*angstrom^2'), symmetry=1, barrier=(17.2878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668697,0.0770865,-0.000100836,7.27333e-08,-2.1164e-11,80377,31.3468], Tmin=(100,'K'), Tmax=(838.384,'K')), NASAPolynomial(coeffs=[11.0489,0.0275637,-1.22349e-05,2.28219e-09,-1.5679e-13,78636.4,-16.902], Tmin=(838.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C=CO[O](22566)',
    structure = SMILES('C=C=C[CH]O[O]'),
    E0 = (334.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.25889,'amu*angstrom^2'), symmetry=1, barrier=(28.9444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25982,'amu*angstrom^2'), symmetry=1, barrier=(28.9658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79467,0.0421231,-2.85308e-05,5.68913e-09,1.05197e-12,40300.5,21.8157], Tmin=(100,'K'), Tmax=(1069.23,'K')), NASAPolynomial(coeffs=[11.7823,0.0154851,-6.20818e-06,1.15293e-09,-8.10381e-14,37551.5,-29.9047], Tmin=(1069.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CCJO)"""),
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
    label = '[CH2][C]=CC=[C][CH2](16912)',
    structure = SMILES('[CH2][C]=CC=[C][CH2]'),
    E0 = (733.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,383.297],'cm^-1')),
        HinderedRotor(inertia=(2.93244,'amu*angstrom^2'), symmetry=1, barrier=(67.4226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191132,'amu*angstrom^2'), symmetry=1, barrier=(67.4298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188775,'amu*angstrom^2'), symmetry=1, barrier=(67.4008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3375.74,'J/mol'), sigma=(5.79513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.28 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84928,0.0398838,-8.2793e-06,-1.72322e-08,9.68702e-12,88323.6,22.6963], Tmin=(100,'K'), Tmax=(960.212,'K')), NASAPolynomial(coeffs=[10.3621,0.0219425,-7.62276e-06,1.3151e-09,-8.95683e-14,85881.1,-22.2335], Tmin=(960.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH][C]([C]=C)O[O](22567)',
    structure = SMILES('[CH2][C]=CC(=[C][CH2])O[O]'),
    E0 = (796.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,2784.69],'cm^-1')),
        HinderedRotor(inertia=(0.199467,'amu*angstrom^2'), symmetry=1, barrier=(4.58613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09504,'amu*angstrom^2'), symmetry=1, barrier=(71.1612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.0962,'amu*angstrom^2'), symmetry=1, barrier=(71.1877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09328,'amu*angstrom^2'), symmetry=1, barrier=(71.1206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823296,0.073808,-9.83795e-05,7.61632e-08,-2.38338e-11,95880,31.119], Tmin=(100,'K'), Tmax=(853.965,'K')), NASAPolynomial(coeffs=[9.04653,0.030172,-1.27425e-05,2.29052e-09,-1.52865e-13,94662.1,-6.16248], Tmin=(853.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(796.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C([CH][C]=C)O[O](22568)',
    structure = SMILES('[CH]=[C]C([CH][C]=C)O[O]'),
    E0 = (907.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,430.011,853.775],'cm^-1')),
        HinderedRotor(inertia=(0.858236,'amu*angstrom^2'), symmetry=1, barrier=(19.7325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858232,'amu*angstrom^2'), symmetry=1, barrier=(19.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0683757,'amu*angstrom^2'), symmetry=1, barrier=(35.3701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858228,'amu*angstrom^2'), symmetry=1, barrier=(19.7324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541791,0.070023,-6.73894e-05,3.22569e-08,-6.08676e-12,109256,34.3935], Tmin=(100,'K'), Tmax=(1284.42,'K')), NASAPolynomial(coeffs=[16.8155,0.0193421,-8.2014e-06,1.53556e-09,-1.0708e-13,105075,-48.1911], Tmin=(1284.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C([C]=C)O[O](22569)',
    structure = SMILES('[CH]=[C][CH]C([C]=C)O[O]'),
    E0 = (907.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,430.011,853.775],'cm^-1')),
        HinderedRotor(inertia=(0.858236,'amu*angstrom^2'), symmetry=1, barrier=(19.7325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858232,'amu*angstrom^2'), symmetry=1, barrier=(19.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0683757,'amu*angstrom^2'), symmetry=1, barrier=(35.3701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858228,'amu*angstrom^2'), symmetry=1, barrier=(19.7324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541791,0.070023,-6.73894e-05,3.22569e-08,-6.08676e-12,109256,34.3935], Tmin=(100,'K'), Tmax=(1284.42,'K')), NASAPolynomial(coeffs=[16.8155,0.0193421,-8.2014e-06,1.53556e-09,-1.0708e-13,105075,-48.1911], Tmin=(1284.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C[C]([C]=C)O[O](22302)',
    structure = SMILES('[CH2][C]=C(C[C]=C)O[O]'),
    E0 = (741.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,299.388,299.397],'cm^-1')),
        HinderedRotor(inertia=(0.0989181,'amu*angstrom^2'), symmetry=1, barrier=(6.28924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0988819,'amu*angstrom^2'), symmetry=1, barrier=(6.28862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0988977,'amu*angstrom^2'), symmetry=1, barrier=(6.28873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0988523,'amu*angstrom^2'), symmetry=1, barrier=(6.28829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949968,0.0745288,-8.67935e-05,4.70221e-08,-3.9734e-12,89326.5,31.7177], Tmin=(100,'K'), Tmax=(601.532,'K')), NASAPolynomial(coeffs=[9.06512,0.0332359,-1.54195e-05,2.93598e-09,-2.03951e-13,88121,-5.21384], Tmin=(601.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH][C](C=C)O[O](22570)',
    structure = SMILES('[CH2][C]=CC(=C[CH2])O[O]'),
    E0 = (558.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955892,0.0679762,-6.89785e-05,4.03894e-08,-9.81854e-12,67272.4,30.1731], Tmin=(100,'K'), Tmax=(985.947,'K')), NASAPolynomial(coeffs=[10.0874,0.0309289,-1.26146e-05,2.27712e-09,-1.54504e-13,65471.8,-13.7517], Tmin=(985.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH][C]=C)O[O](22571)',
    structure = SMILES('[CH]=CC([CH][C]=C)O[O]'),
    E0 = (669.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,527.971,772.109],'cm^-1')),
        HinderedRotor(inertia=(0.0548253,'amu*angstrom^2'), symmetry=1, barrier=(23.078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00421,'amu*angstrom^2'), symmetry=1, barrier=(23.0887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00483,'amu*angstrom^2'), symmetry=1, barrier=(23.1029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00394,'amu*angstrom^2'), symmetry=1, barrier=(23.0824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494692,0.0663909,-4.61368e-05,7.86406e-09,2.64697e-12,80655.8,34.0877], Tmin=(100,'K'), Tmax=(1045.7,'K')), NASAPolynomial(coeffs=[17.5922,0.0204955,-8.28184e-06,1.56804e-09,-1.12325e-13,76013.6,-54.2611], Tmin=(1045.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]C([C]=C)O[O](22572)',
    structure = SMILES('[CH]=C[CH]C([C]=C)O[O]'),
    E0 = (669.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,527.971,772.109],'cm^-1')),
        HinderedRotor(inertia=(0.0548253,'amu*angstrom^2'), symmetry=1, barrier=(23.078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00421,'amu*angstrom^2'), symmetry=1, barrier=(23.0887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00483,'amu*angstrom^2'), symmetry=1, barrier=(23.1029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00394,'amu*angstrom^2'), symmetry=1, barrier=(23.0824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494692,0.0663909,-4.61368e-05,7.86406e-09,2.64697e-12,80655.8,34.0877], Tmin=(100,'K'), Tmax=(1045.7,'K')), NASAPolynomial(coeffs=[17.5922,0.0204955,-8.28184e-06,1.56804e-09,-1.12325e-13,76013.6,-54.2611], Tmin=(1045.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH][C]([C]=C)OO(22573)',
    structure = SMILES('[CH2][C]=CC(=[C][CH2])OO'),
    E0 = (644.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707475,0.0755257,-8.65014e-05,5.66243e-08,-1.52776e-11,77603.3,30.7844], Tmin=(100,'K'), Tmax=(894.57,'K')), NASAPolynomial(coeffs=[10.3756,0.0322937,-1.40081e-05,2.59758e-09,-1.78586e-13,75873.6,-14.7814], Tmin=(894.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][C]([CH]C=C)O[O](22574)',
    structure = SMILES('[CH2][C]=C(C=C[CH2])O[O]'),
    E0 = (558.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955892,0.0679762,-6.89785e-05,4.03894e-08,-9.81854e-12,67272.4,30.1731], Tmin=(100,'K'), Tmax=(985.947,'K')), NASAPolynomial(coeffs=[10.0874,0.0309289,-1.26146e-05,2.27712e-09,-1.54504e-13,65471.8,-13.7517], Tmin=(985.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556652,0.080901,-0.000105444,7.68688e-08,-2.28706e-11,100809,33.5511], Tmin=(100,'K'), Tmax=(816.381,'K')), NASAPolynomial(coeffs=[10.6219,0.0315779,-1.48073e-05,2.84365e-09,-1.98929e-13,99165.9,-12.9648], Tmin=(816.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(837.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(C[C]=C)O[O](22303)',
    structure = SMILES('[CH]=[C]C(C[C]=C)O[O]'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556652,0.080901,-0.000105444,7.68688e-08,-2.28706e-11,100809,33.5511], Tmin=(100,'K'), Tmax=(816.381,'K')), NASAPolynomial(coeffs=[10.6219,0.0315779,-1.48073e-05,2.84365e-09,-1.98929e-13,99165.9,-12.9648], Tmin=(816.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(837.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH][C]=C)OO(22575)',
    structure = SMILES('[CH]=[C]C([CH][C]=C)OO'),
    E0 = (755.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,474.187,474.552],'cm^-1')),
        HinderedRotor(inertia=(0.171516,'amu*angstrom^2'), symmetry=1, barrier=(27.3239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171144,'amu*angstrom^2'), symmetry=1, barrier=(27.3229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170648,'amu*angstrom^2'), symmetry=1, barrier=(27.3228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.586668,'amu*angstrom^2'), symmetry=1, barrier=(27.3237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170924,'amu*angstrom^2'), symmetry=1, barrier=(27.3216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0787632,0.0758282,-6.97855e-05,3.12353e-08,-5.46564e-12,90994,35.305], Tmin=(100,'K'), Tmax=(1387.89,'K')), NASAPolynomial(coeffs=[19.9836,0.0184609,-7.78414e-06,1.45326e-09,-1.01004e-13,85468.9,-67.2489], Tmin=(1387.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH]C=C)O[O](22576)',
    structure = SMILES('[CH]=[C]C([CH]C=C)O[O]'),
    E0 = (669.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,527.971,772.109],'cm^-1')),
        HinderedRotor(inertia=(0.0548253,'amu*angstrom^2'), symmetry=1, barrier=(23.078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00421,'amu*angstrom^2'), symmetry=1, barrier=(23.0887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00483,'amu*angstrom^2'), symmetry=1, barrier=(23.1029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00394,'amu*angstrom^2'), symmetry=1, barrier=(23.0824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494692,0.0663909,-4.61368e-05,7.86406e-09,2.64697e-12,80655.8,34.0877], Tmin=(100,'K'), Tmax=(1045.7,'K')), NASAPolynomial(coeffs=[17.5922,0.0204955,-8.28184e-06,1.56804e-09,-1.12325e-13,76013.6,-54.2611], Tmin=(1045.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(C=C)O[O](22577)',
    structure = SMILES('[CH]=[C][CH]C(C=C)O[O]'),
    E0 = (669.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,771.085],'cm^-1')),
        HinderedRotor(inertia=(1.00418,'amu*angstrom^2'), symmetry=1, barrier=(23.0881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00422,'amu*angstrom^2'), symmetry=1, barrier=(23.0891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11672,'amu*angstrom^2'), symmetry=1, barrier=(23.0874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0547007,'amu*angstrom^2'), symmetry=1, barrier=(23.0875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494692,0.0663909,-4.61368e-05,7.86406e-09,2.64697e-12,80655.8,34.0877], Tmin=(100,'K'), Tmax=(1045.7,'K')), NASAPolynomial(coeffs=[17.5922,0.0204955,-8.28184e-06,1.56804e-09,-1.12325e-13,76013.6,-54.2611], Tmin=(1045.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C([C]=C)OO(22578)',
    structure = SMILES('[CH]=[C][CH]C([C]=C)OO'),
    E0 = (755.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,474.187,474.552],'cm^-1')),
        HinderedRotor(inertia=(0.171516,'amu*angstrom^2'), symmetry=1, barrier=(27.3239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171144,'amu*angstrom^2'), symmetry=1, barrier=(27.3229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170648,'amu*angstrom^2'), symmetry=1, barrier=(27.3228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.586668,'amu*angstrom^2'), symmetry=1, barrier=(27.3237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170924,'amu*angstrom^2'), symmetry=1, barrier=(27.3216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0787632,0.0758282,-6.97855e-05,3.12353e-08,-5.46564e-12,90994,35.305], Tmin=(100,'K'), Tmax=(1387.89,'K')), NASAPolynomial(coeffs=[19.9836,0.0184609,-7.78414e-06,1.45326e-09,-1.01004e-13,85468.9,-67.2489], Tmin=(1387.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJC(O)C=C) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH2]C1=CC(O[O])C1=C(22579)',
    structure = SMILES('C=C1[CH]C(O[O])C1=C'),
    E0 = (307.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18035,0.0404599,4.03516e-05,-8.86777e-08,3.84683e-11,37045.2,24.9829], Tmin=(100,'K'), Tmax=(958.834,'K')), NASAPolynomial(coeffs=[19.9117,0.0158315,-4.8364e-06,9.48391e-10,-7.68105e-14,30993.3,-77.4251], Tmin=(958.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH2]C1=CC([C]=C)OO1(22580)',
    structure = SMILES('C=[C]C1[CH]C(=C)OO1'),
    E0 = (354.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76939,0.0210978,9.66878e-05,-1.45933e-07,5.81086e-11,42801.1,28.3408], Tmin=(100,'K'), Tmax=(959.964,'K')), NASAPolynomial(coeffs=[20.0459,0.0154635,-4.70086e-06,1.00369e-09,-8.68868e-14,36042.7,-76.0102], Tmin=(959.964,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC([C]=C)O1(22581)',
    structure = SMILES('C=[C]C1[CH]C(=C)O1'),
    E0 = (319.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59616,0.0222267,9.95011e-05,-1.65067e-07,7.09683e-11,38601.1,21.0371], Tmin=(100,'K'), Tmax=(914.303,'K')), NASAPolynomial(coeffs=[24.9776,-0.000307178,5.61967e-06,-1.20263e-09,7.44461e-14,30991.9,-107.901], Tmin=(914.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CCJC(O)C=C) + radical(Cds_S)"""),
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
    label = '[CH2][C]=CC#C[CH2](17826)',
    structure = SMILES('[CH2][C]=CC#C[CH2]'),
    E0 = (708.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2100,2250,500,550,322.101],'cm^-1')),
        HinderedRotor(inertia=(0.00162484,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664835,'amu*angstrom^2'), symmetry=1, barrier=(48.946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664825,'amu*angstrom^2'), symmetry=1, barrier=(48.9461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87382,0.0430323,-3.17318e-05,1.24373e-08,-2.01431e-12,85252.2,20.9301], Tmin=(100,'K'), Tmax=(1434.79,'K')), NASAPolynomial(coeffs=[10.0068,0.0203588,-8.02804e-06,1.42354e-09,-9.52669e-14,82918.4,-21.243], Tmin=(1434.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CTCC=CCJ) + radical(Propargyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC(=C=C)O[O](22582)',
    structure = SMILES('[CH2]C=CC(=C=C)O[O]'),
    E0 = (379.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788399,0.0714444,-7.48016e-05,4.32748e-08,-1.02254e-11,45714.6,27.7965], Tmin=(100,'K'), Tmax=(1019.17,'K')), NASAPolynomial(coeffs=[11.6574,0.0287864,-1.20184e-05,2.20673e-09,-1.51533e-13,43499.1,-24.8465], Tmin=(1019.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C#CC([C]=C)O[O](22583)',
    structure = SMILES('[CH2]C#CC([C]=C)O[O]'),
    E0 = (652.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,2100,2250,500,550,269.397,269.412],'cm^-1')),
        HinderedRotor(inertia=(0.218226,'amu*angstrom^2'), symmetry=1, barrier=(11.2404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00232275,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09028,'amu*angstrom^2'), symmetry=1, barrier=(107.659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09055,'amu*angstrom^2'), symmetry=1, barrier=(107.659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925318,0.0707961,-8.46902e-05,5.69517e-08,-1.56519e-11,78578.8,31.0349], Tmin=(100,'K'), Tmax=(881.557,'K')), NASAPolynomial(coeffs=[10.2813,0.028344,-1.24563e-05,2.32554e-09,-1.60446e-13,76929.3,-12.9227], Tmin=(881.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(Propargyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]C([C]=C)O[O](22584)',
    structure = SMILES('[CH2]C=[C]C([C]=C)O[O]'),
    E0 = (725.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,334.331,334.516],'cm^-1')),
        HinderedRotor(inertia=(0.121483,'amu*angstrom^2'), symmetry=1, barrier=(9.64996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121642,'amu*angstrom^2'), symmetry=1, barrier=(9.65082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121582,'amu*angstrom^2'), symmetry=1, barrier=(9.65155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369598,'amu*angstrom^2'), symmetry=1, barrier=(29.3496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771081,0.0762149,-9.13859e-05,6.23742e-08,-1.76797e-11,87379.2,32.0619], Tmin=(100,'K'), Tmax=(848.767,'K')), NASAPolynomial(coeffs=[9.83616,0.0334968,-1.58971e-05,3.0857e-09,-2.17844e-13,85840.3,-10.1861], Tmin=(848.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(725.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(C=C)O[O](22585)',
    structure = SMILES('[CH2][C]=[C]C(C=C)O[O]'),
    E0 = (725.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,334.479,334.695],'cm^-1')),
        HinderedRotor(inertia=(0.121657,'amu*angstrom^2'), symmetry=1, barrier=(9.65045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121558,'amu*angstrom^2'), symmetry=1, barrier=(9.65135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121675,'amu*angstrom^2'), symmetry=1, barrier=(9.65053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369764,'amu*angstrom^2'), symmetry=1, barrier=(29.3496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771081,0.0762149,-9.13859e-05,6.23742e-08,-1.76797e-11,87379.2,32.0619], Tmin=(100,'K'), Tmax=(848.767,'K')), NASAPolynomial(coeffs=[9.83616,0.0334968,-1.58971e-05,3.0857e-09,-2.17844e-13,85840.3,-10.1861], Tmin=(848.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(725.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C([C]=[C]C)O[O](22586)',
    structure = SMILES('C=[C]C([C]=[C]C)O[O]'),
    E0 = (811.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.368527,'amu*angstrom^2'), symmetry=1, barrier=(8.47315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368755,'amu*angstrom^2'), symmetry=1, barrier=(8.4784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36863,'amu*angstrom^2'), symmetry=1, barrier=(8.47552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368685,'amu*angstrom^2'), symmetry=1, barrier=(8.4768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.485193,0.0865855,-0.000137687,1.24259e-07,-4.4107e-11,97769.8,33.2628], Tmin=(100,'K'), Tmax=(831.421,'K')), NASAPolynomial(coeffs=[6.98464,0.0379245,-1.85186e-05,3.54489e-09,-2.44229e-13,97290.1,6.72162], Tmin=(831.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(811.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([C]=C)OO(22587)',
    structure = SMILES('[CH2][C]=[C]C([C]=C)OO'),
    E0 = (811.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,267.429,269.796],'cm^-1')),
        HinderedRotor(inertia=(0.272954,'amu*angstrom^2'), symmetry=1, barrier=(14.229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00233223,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275354,'amu*angstrom^2'), symmetry=1, barrier=(14.2212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695059,'amu*angstrom^2'), symmetry=1, barrier=(35.9061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698817,'amu*angstrom^2'), symmetry=1, barrier=(35.939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434379,0.084896,-0.000113313,8.5027e-08,-2.62291e-11,97713.8,32.984], Tmin=(100,'K'), Tmax=(785.54,'K')), NASAPolynomial(coeffs=[10.3331,0.0344932,-1.70726e-05,3.35368e-09,-2.37514e-13,96158.6,-12.3829], Tmin=(785.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(811.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][C](C=[C]C)O[O](22588)',
    structure = SMILES('[CH2][C]=C(C=[C]C)O[O]'),
    E0 = (678.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.29319,'amu*angstrom^2'), symmetry=1, barrier=(6.74102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292887,'amu*angstrom^2'), symmetry=1, barrier=(6.73404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292879,'amu*angstrom^2'), symmetry=1, barrier=(6.73387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901973,'amu*angstrom^2'), symmetry=1, barrier=(20.7381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549574,0.0826638,-0.000121653,1.0251e-07,-3.42384e-11,81688.3,30.5423], Tmin=(100,'K'), Tmax=(856.029,'K')), NASAPolynomial(coeffs=[8.23567,0.0343691,-1.53348e-05,2.81711e-09,-1.89587e-13,80826,-2.69459], Tmin=(856.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(678.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(C=[C]C)O[O](22589)',
    structure = SMILES('[CH]=[C]C(C=[C]C)O[O]'),
    E0 = (821.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.450396,'amu*angstrom^2'), symmetry=1, barrier=(10.3555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449588,'amu*angstrom^2'), symmetry=1, barrier=(10.3369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449783,'amu*angstrom^2'), symmetry=1, barrier=(10.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449579,'amu*angstrom^2'), symmetry=1, barrier=(10.3367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461086,0.0852917,-0.000128701,1.11124e-07,-3.84682e-11,98885.4,33.2801], Tmin=(100,'K'), Tmax=(810.31,'K')), NASAPolynomial(coeffs=[8.36577,0.0356812,-1.72616e-05,3.31054e-09,-2.2918e-13,97952.1,-1.04754], Tmin=(810.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(821.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (660.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (660.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (910.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1107.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1202.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1175.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (663.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (668.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (668.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (791.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (791.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (688.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (683.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (723.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (825.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (860.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (894.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (934.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (878.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (660.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (725.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1008.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1119.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1119.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (900.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (823.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (774.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (824.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (825.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (785.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (982.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1046.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (788.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (978.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (978.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (906.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (754.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (668.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (667.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (743.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (838.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (738.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (879.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (798.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (929.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (901.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (973.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (895.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (806.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (953.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C3H3(5450)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['O2(S)(666)', '[CH2][C]=CC=C=C(16910)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', 'C=[C][CH]C([O])[C]=C(22553)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C][CH2](16918)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H2CC(T)(1341)', 'C=[C][CH][CH]O[O](22554)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=[C][C]C([C]=C)O[O](22555)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C]C1C(=C)C1O[O](22556)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C]C1OOC1[C]=C(22557)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C][CH]C1OOC1=C(22558)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['O(4)', 'C=[C]C1OC1[C]=C(22559)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['O(4)', 'C=[C][CH]C1OC1=C(22560)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO_S;Y_rad_intra;OOJ] for rate rule [R2OO_S;Cd_rad_out;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C]CC(=C=C)O[O](22300)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C]C=C(C=C)O[O](22561)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C]C=C([C]=C)OO(22562)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C=[C]C=C([C]=C)O[O](22563)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C#CC([CH][C]=C)O[O](22564)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', 'C#C[CH]C([C]=C)O[O](22565)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H2CC(T)(1341)', 'C=[C]C=CO[O](22566)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][CH2](16918)', 'C=C=CO[O](16806)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00817,'m^3/(mol*s)'), n=1.99965, Ea=(13.9682,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;YJ] for rate rule [Cds_Ca;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O2(2)', '[CH2][C]=CC=C=C(16910)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.000334811,'m^3/(mol*s)'), n=2.98833, Ea=(125.168,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O2b]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from 122.5 to 125.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['O2(2)', '[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.4339e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'C=[C][CH][C]([C]=C)O[O](22567)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH]=[C]C([CH][C]=C)O[O](22568)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH]=[C][CH]C([C]=C)O[O](22569)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C[C]([C]=C)O[O](22302)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C][CH][C](C=C)O[O](22570)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC([CH][C]=C)O[O](22571)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C[CH]C([C]=C)O[O](22572)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C][CH][C]([C]=C)OO(22573)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.83109e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['C=[C][C]([CH]C=C)O[O](22574)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.30234e+06,'s^-1'), n=1.68744, Ea=(125.264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]CC([C]=C)O[O](21285)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=[C]C(C[C]=C)O[O](22303)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=[C]C([CH][C]=C)OO(22575)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]C([CH]C=C)O[O](22576)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C][CH]C(C=C)O[O](22577)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=[C][CH]C([C]=C)OO(22578)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['[CH2][C]=CC(=C)[CH]O[O](20233)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['[CH2]C1=CC(O[O])C1=C(22579)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['[CH2]C1=CC([C]=C)OO1(22580)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['O(4)', '[CH2]C1=CC([C]=C)O1(22581)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO] for rate rule [R3OO_DS;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['HO2(9)', '[CH2][C]=CC#C[CH2](17826)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    products = ['[CH2]C=CC(=C=C)O[O](22582)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)', '[CH2]C#CC([C]=C)O[O](22583)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C3H3(5450)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=[C]C([C]=C)O[O](22584)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][C]=[C]C(C=C)O[O](22585)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=[C]C([C]=[C]C)O[O](22586)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][C]=[C]C([C]=C)OO(22587)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;O_H_out] for rate rule [R4H_SSS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=[C][C](C=[C]C)O[O](22588)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(190215,'s^-1'), n=1.95446, Ea=(128.468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=[C]C(C=[C]C)O[O](22589)'],
    products = ['[CH2][C]=CC([C]=C)O[O](21296)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #5018',
    isomers = [
        '[CH2][C]=CC([C]=C)O[O](21296)',
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
    label = 'PDepNetwork #5018',
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


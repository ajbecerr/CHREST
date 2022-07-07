species(
    label = '[CH]=[C]CC=[C]CO[O](19500)',
    structure = SMILES('[CH]=[C]CC=[C]CO[O]'),
    E0 = (845.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.307449,'amu*angstrom^2'), symmetry=1, barrier=(7.06887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307457,'amu*angstrom^2'), symmetry=1, barrier=(7.06903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307414,'amu*angstrom^2'), symmetry=1, barrier=(7.06805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307429,'amu*angstrom^2'), symmetry=1, barrier=(7.06841,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563494,0.0834111,-0.000126883,1.10744e-07,-3.83322e-11,101852,33.5899], Tmin=(100,'K'), Tmax=(835.565,'K')), NASAPolynomial(coeffs=[7.56703,0.0360318,-1.69609e-05,3.20054e-09,-2.1897e-13,101165,3.95365], Tmin=(835.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = 'C#CCO[O](16808)',
    structure = SMILES('C#CCO[O]'),
    E0 = (242.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.595127,'amu*angstrom^2'), symmetry=1, barrier=(13.6831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74124,'amu*angstrom^2'), symmetry=1, barrier=(63.0265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3504.08,'J/mol'), sigma=(5.7666,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=547.33 K, Pc=41.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29019,0.0394598,-5.33344e-05,4.08815e-08,-1.22872e-11,29219.5,16.1946], Tmin=(100,'K'), Tmax=(942.252,'K')), NASAPolynomial(coeffs=[6.8169,0.0148149,-5.46001e-06,8.95137e-10,-5.57485e-14,28607.4,-4.09622], Tmin=(942.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C([CH2])C=[C]CO[O](19498)',
    structure = SMILES('[CH]C(=C)C=[C]CO[O]'),
    E0 = (699.205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445273,0.081847,-9.32577e-05,6.15055e-08,-1.684e-11,84219.7,30.0304], Tmin=(100,'K'), Tmax=(879.583,'K')), NASAPolynomial(coeffs=[10.3995,0.0365758,-1.60491e-05,2.98248e-09,-2.05115e-13,82468.7,-16.7153], Tmin=(879.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=[C]CC=C=C(16903)',
    structure = SMILES('[CH]=[C]CC=C=C'),
    E0 = (714.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.793161,'amu*angstrom^2'), symmetry=1, barrier=(18.2363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.795875,'amu*angstrom^2'), symmetry=1, barrier=(18.2987,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68641,0.0494647,-4.02886e-05,1.7498e-08,-3.14662e-12,85990.8,22.8642], Tmin=(100,'K'), Tmax=(1298.36,'K')), NASAPolynomial(coeffs=[10.2755,0.0230035,-9.71798e-06,1.80111e-09,-1.24189e-13,83760.4,-20.816], Tmin=(1298.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC=[C]C[O](20309)',
    structure = SMILES('[CH]=[C]CC=[C]C[O]'),
    E0 = (848.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,196.267,196.3,2282.2],'cm^-1')),
        HinderedRotor(inertia=(0.216447,'amu*angstrom^2'), symmetry=1, barrier=(5.91647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216448,'amu*angstrom^2'), symmetry=1, barrier=(5.9168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216381,'amu*angstrom^2'), symmetry=1, barrier=(5.91688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33702,0.0673553,-0.000103533,9.86353e-08,-3.70126e-11,102087,30.1519], Tmin=(100,'K'), Tmax=(830.793,'K')), NASAPolynomial(coeffs=[3.27385,0.039281,-1.89936e-05,3.63264e-09,-2.50576e-13,102412,25.061], Tmin=(830.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(848.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH2]C=[C]CO[O](19545)',
    structure = SMILES('[CH2]C=[C]CO[O]'),
    E0 = (429.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.175324,'amu*angstrom^2'), symmetry=1, barrier=(4.03105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174292,'amu*angstrom^2'), symmetry=1, barrier=(4.00731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942554,'amu*angstrom^2'), symmetry=1, barrier=(21.6712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76322,0.0527063,-6.55699e-05,5.10847e-08,-1.68387e-11,51771.8,23.0334], Tmin=(100,'K'), Tmax=(766.437,'K')), NASAPolynomial(coeffs=[6.45998,0.0268127,-1.21894e-05,2.30126e-09,-1.59132e-13,51092.4,1.88831], Tmin=(766.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[C]=CC[C]=[CH](17531)',
    structure = SMILES('[C]=CC[C]=[CH]'),
    E0 = (1131.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,239.393],'cm^-1')),
        HinderedRotor(inertia=(0.277767,'amu*angstrom^2'), symmetry=1, barrier=(10.7298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257782,'amu*angstrom^2'), symmetry=1, barrier=(10.7284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30755,0.0399092,-4.49451e-05,2.98959e-08,-8.40636e-12,136181,20.1622], Tmin=(100,'K'), Tmax=(849.57,'K')), NASAPolynomial(coeffs=[6.62319,0.0195903,-9.07074e-06,1.74536e-09,-1.22727e-13,135448,0.0452128], Tmin=(849.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1131.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
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
    label = '[C]=[C]CC=[C]CO[O](20310)',
    structure = SMILES('[C]=[C]CC=[C]CO[O]'),
    E0 = (1156.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,2270.13],'cm^-1')),
        HinderedRotor(inertia=(0.416753,'amu*angstrom^2'), symmetry=1, barrier=(9.58196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416704,'amu*angstrom^2'), symmetry=1, barrier=(9.58085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416929,'amu*angstrom^2'), symmetry=1, barrier=(9.58603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416485,'amu*angstrom^2'), symmetry=1, barrier=(9.57581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576034,0.086572,-0.000148569,1.38391e-07,-4.9326e-11,139253,33.5185], Tmin=(100,'K'), Tmax=(856.992,'K')), NASAPolynomial(coeffs=[6.10925,0.036301,-1.77935e-05,3.37434e-09,-2.29732e-13,139203,12.9151], Tmin=(856.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1156.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC=C1COO1(20311)',
    structure = SMILES('[CH]=[C]CC=C1COO1'),
    E0 = (622.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21047,0.050887,-1.64918e-05,-1.15404e-08,6.69497e-12,74967.1,30.5366], Tmin=(100,'K'), Tmax=(1110,'K')), NASAPolynomial(coeffs=[13.5287,0.026975,-1.18511e-05,2.29271e-09,-1.63889e-13,70970.9,-35.8601], Tmin=(1110,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC=C1CO[O](20312)',
    structure = SMILES('[CH]=C1CC=C1CO[O]'),
    E0 = (499.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.718,0.0614318,-3.42744e-05,-5.2026e-09,7.98809e-12,60216.2,26.572], Tmin=(100,'K'), Tmax=(977.587,'K')), NASAPolynomial(coeffs=[16.8402,0.0195002,-6.81502e-06,1.22182e-09,-8.6582e-14,55915.5,-56.7171], Tmin=(977.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC=C1CO1(20313)',
    structure = SMILES('[CH]=[C]CC=C1CO1'),
    E0 = (555.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13814,0.0497454,-1.00562e-05,-3.11226e-08,1.79086e-11,66965.2,24.0848], Tmin=(100,'K'), Tmax=(947.382,'K')), NASAPolynomial(coeffs=[17.5414,0.0130698,-3.57429e-06,6.17496e-10,-4.65225e-14,62395,-61.8819], Tmin=(947.382,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC=CCO[O](20314)',
    structure = SMILES('C#CC=C[CH]CO[O]'),
    E0 = (392.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.401029,0.0783247,-8.04226e-05,4.27059e-08,-9.10409e-12,47353.4,25.9523], Tmin=(100,'K'), Tmax=(1130.66,'K')), NASAPolynomial(coeffs=[15.0226,0.0265975,-1.17987e-05,2.24368e-09,-1.57554e-13,44047,-46.3839], Tmin=(1130.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=[C]CC=C=COO(20315)',
    structure = SMILES('[CH]=[C]CC=C=COO'),
    E0 = (629.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449854,0.0793754,-8.95335e-05,5.27847e-08,-1.24655e-11,75894.6,31.2125], Tmin=(100,'K'), Tmax=(1026.97,'K')), NASAPolynomial(coeffs=[14.1422,0.0260443,-1.16379e-05,2.21817e-09,-1.55841e-13,73082.2,-35.2098], Tmin=(1026.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC=C=CO[O](19507)',
    structure = SMILES('[CH]=CCC=C=CO[O]'),
    E0 = (544.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.614046,0.0728255,-7.55217e-05,4.11089e-08,-8.95118e-12,65567.3,30.9027], Tmin=(100,'K'), Tmax=(1112.64,'K')), NASAPolynomial(coeffs=[14.1374,0.0242077,-9.97683e-06,1.83532e-09,-1.26635e-13,62558,-35.7827], Tmin=(1112.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C#C[CH]C=[C]CO[O](20280)',
    structure = SMILES('[CH]=C=CC=[C]CO[O]'),
    E0 = (661.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,540,610,2055,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.792744,'amu*angstrom^2'), symmetry=1, barrier=(18.2267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796024,'amu*angstrom^2'), symmetry=1, barrier=(18.3022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793139,'amu*angstrom^2'), symmetry=1, barrier=(18.2358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438553,0.0813041,-0.000112382,8.38726e-08,-2.47668e-11,79645.4,29.8211], Tmin=(100,'K'), Tmax=(870.475,'K')), NASAPolynomial(coeffs=[12.1376,0.0247253,-1.00268e-05,1.76127e-09,-1.15748e-13,77715.5,-24.3835], Tmin=(870.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]CC=C=CO[O](20316)',
    structure = SMILES('[CH]=[C]CC=C=CO[O]'),
    E0 = (781.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.647565,'amu*angstrom^2'), symmetry=1, barrier=(14.8888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647624,'amu*angstrom^2'), symmetry=1, barrier=(14.8901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646982,'amu*angstrom^2'), symmetry=1, barrier=(14.8754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650742,0.0765139,-9.67577e-05,6.53742e-08,-1.76586e-11,94167.8,31.2513], Tmin=(100,'K'), Tmax=(903.971,'K')), NASAPolynomial(coeffs=[12.3475,0.0247556,-1.08705e-05,2.03197e-09,-1.40346e-13,92053.1,-23.9974], Tmin=(903.971,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(781.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC#CCO[O](20317)',
    structure = SMILES('[CH]=[C]CC#CCO[O]'),
    E0 = (771.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2100,2250,500,550,1685,370,3120,650,792.5,1650,269.61,269.61],'cm^-1')),
        HinderedRotor(inertia=(0.130106,'amu*angstrom^2'), symmetry=1, barrier=(6.71113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130105,'amu*angstrom^2'), symmetry=1, barrier=(6.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130105,'amu*angstrom^2'), symmetry=1, barrier=(6.71112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21997,'amu*angstrom^2'), symmetry=1, barrier=(62.9289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.836021,0.0746319,-0.000105603,8.55183e-08,-2.75162e-11,92870.3,30.9207], Tmin=(100,'K'), Tmax=(872.496,'K')), NASAPolynomial(coeffs=[8.5491,0.0299356,-1.27113e-05,2.27757e-09,-1.51026e-13,91879.7,-3.20217], Tmin=(872.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CO[O](16805)',
    structure = SMILES('[CH]=[C]CO[O]'),
    E0 = (567.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.348758,'amu*angstrom^2'), symmetry=1, barrier=(8.01863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349832,'amu*angstrom^2'), symmetry=1, barrier=(8.04334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92348,0.0566508,-0.000115722,1.14457e-07,-4.10713e-11,68272.1,19.3685], Tmin=(100,'K'), Tmax=(892.364,'K')), NASAPolynomial(coeffs=[4.05419,0.0203884,-9.86661e-06,1.83106e-09,-1.20972e-13,68955.4,15.2908], Tmin=(892.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC=[C][CH2](16907)',
    structure = SMILES('[CH]=[C]CC=[C][CH2]'),
    E0 = (927.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0544244,'amu*angstrom^2'), symmetry=1, barrier=(10.4638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135165,'amu*angstrom^2'), symmetry=1, barrier=(25.0247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348043,'amu*angstrom^2'), symmetry=1, barrier=(10.5051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3289.69,'J/mol'), sigma=(5.72141,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.84 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68217,0.0511066,-4.45835e-05,2.12564e-08,-4.22835e-12,111576,24.5816], Tmin=(100,'K'), Tmax=(1179.05,'K')), NASAPolynomial(coeffs=[9.59155,0.0242736,-1.04464e-05,1.95443e-09,-1.35673e-13,109711,-14.8793], Tmin=(1179.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(927.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=CC=[C]CO[O](20256)',
    structure = SMILES('[CH]=C=C[CH][C]CO[O]'),
    E0 = (907.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.79104,'amu*angstrom^2'), symmetry=1, barrier=(64.1714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63723,'amu*angstrom^2'), symmetry=1, barrier=(14.6512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79163,'amu*angstrom^2'), symmetry=1, barrier=(64.185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0142611,'amu*angstrom^2'), symmetry=1, barrier=(7.05605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444679,0.081834,-0.000110952,8.2316e-08,-2.44478e-11,109303,29.5973], Tmin=(100,'K'), Tmax=(824.722,'K')), NASAPolynomial(coeffs=[11.6845,0.0273155,-1.17869e-05,2.14996e-09,-1.45162e-13,107450,-22.4617], Tmin=(824.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]CC=[C][CH]O[O](20318)',
    structure = SMILES('[CH]=[C]CC=[C][CH]O[O]'),
    E0 = (963.173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1670,1700,300,440,3120,650,792.5,1650,322.857,322.86],'cm^-1')),
        HinderedRotor(inertia=(0.00161729,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151526,'amu*angstrom^2'), symmetry=1, barrier=(11.2079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151522,'amu*angstrom^2'), symmetry=1, barrier=(11.208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459394,'amu*angstrom^2'), symmetry=1, barrier=(33.9808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940607,0.0701778,-7.8326e-05,4.74054e-08,-1.17181e-11,115951,32.8452], Tmin=(100,'K'), Tmax=(973.85,'K')), NASAPolynomial(coeffs=[11.3262,0.0275196,-1.26196e-05,2.4244e-09,-1.70732e-13,113928,-16.9838], Tmin=(973.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(963.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C[C]=[C]CO[O](20319)',
    structure = SMILES('[CH]=[C]C[C]=[C]CO[O]'),
    E0 = (1083.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1685,1700,300,370,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.287139,'amu*angstrom^2'), symmetry=1, barrier=(6.6019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287278,'amu*angstrom^2'), symmetry=1, barrier=(6.60509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2874,'amu*angstrom^2'), symmetry=1, barrier=(6.60789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287097,'amu*angstrom^2'), symmetry=1, barrier=(6.60093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504748,0.0883257,-0.000152811,1.41508e-07,-4.99193e-11,130456,34.2736], Tmin=(100,'K'), Tmax=(865.261,'K')), NASAPolynomial(coeffs=[6.51466,0.0352915,-1.70968e-05,3.21554e-09,-2.17443e-13,130361,11.6112], Tmin=(865.261,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1083.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]C=[C]CO[O](20281)',
    structure = SMILES('[CH]C=CC=[C]CO[O]'),
    E0 = (700.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650906,0.0773073,-8.28202e-05,5.22105e-08,-1.38745e-11,84400.6,30.4987], Tmin=(100,'K'), Tmax=(898.626,'K')), NASAPolynomial(coeffs=[9.55849,0.0376578,-1.6637e-05,3.11121e-09,-2.15066e-13,82799.7,-11.5232], Tmin=(898.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC=C[CH]O[O](20320)',
    structure = SMILES('[CH]=[C]CC=C[CH]O[O]'),
    E0 = (725.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793782,0.0678201,-6.18318e-05,2.93149e-08,-5.62335e-12,87354.8,32.8886], Tmin=(100,'K'), Tmax=(1243.49,'K')), NASAPolynomial(coeffs=[13.8988,0.0256654,-1.09824e-05,2.05377e-09,-1.42694e-13,84095.6,-33.1916], Tmin=(1243.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(725.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C[C]=CCO[O](20321)',
    structure = SMILES('[CH]=[C]C[C]=CCO[O]'),
    E0 = (845.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.307449,'amu*angstrom^2'), symmetry=1, barrier=(7.06887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307457,'amu*angstrom^2'), symmetry=1, barrier=(7.06903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307414,'amu*angstrom^2'), symmetry=1, barrier=(7.06805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307429,'amu*angstrom^2'), symmetry=1, barrier=(7.06841,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563494,0.0834111,-0.000126883,1.10744e-07,-3.83322e-11,101852,33.5899], Tmin=(100,'K'), Tmax=(835.565,'K')), NASAPolynomial(coeffs=[7.56703,0.0360318,-1.69609e-05,3.20054e-09,-2.1897e-13,101165,3.95365], Tmin=(835.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC=[C][CH]OO(20322)',
    structure = SMILES('[CH]=[C]CC=[C][CH]OO'),
    E0 = (811.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687311,0.0737169,-7.37043e-05,3.84378e-08,-8.14668e-12,97679.6,32.9897], Tmin=(100,'K'), Tmax=(1127.17,'K')), NASAPolynomial(coeffs=[13.5512,0.0280669,-1.29552e-05,2.50782e-09,-1.77653e-13,94779.6,-30.6111], Tmin=(1127.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(811.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC[C]=[C]CO[O](20323)',
    structure = SMILES('[CH]=CC[C]=[C]CO[O]'),
    E0 = (845.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.307449,'amu*angstrom^2'), symmetry=1, barrier=(7.06887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307457,'amu*angstrom^2'), symmetry=1, barrier=(7.06903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307414,'amu*angstrom^2'), symmetry=1, barrier=(7.06805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307429,'amu*angstrom^2'), symmetry=1, barrier=(7.06841,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563494,0.0834111,-0.000126883,1.10744e-07,-3.83322e-11,101852,33.5899], Tmin=(100,'K'), Tmax=(835.565,'K')), NASAPolynomial(coeffs=[7.56703,0.0360318,-1.69609e-05,3.20054e-09,-2.1897e-13,101165,3.95365], Tmin=(835.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C=CCO[O](20282)',
    structure = SMILES('[CH]=C=C[CH][CH]CO[O]'),
    E0 = (659.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520326,0.0808472,-0.000107607,8.23704e-08,-2.56537e-11,79496.5,31.7242], Tmin=(100,'K'), Tmax=(812.522,'K')), NASAPolynomial(coeffs=[10.0065,0.0321552,-1.40386e-05,2.58054e-09,-1.75063e-13,78020.7,-11.6668], Tmin=(812.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCJCOOH) + radical(Allyl_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2][C]=CC=[C]CO[O](19509)',
    structure = SMILES('[CH2][C]=CC=[C]CO[O]'),
    E0 = (685.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.2086,'amu*angstrom^2'), symmetry=1, barrier=(4.79612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799987,'amu*angstrom^2'), symmetry=1, barrier=(18.3933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208521,'amu*angstrom^2'), symmetry=1, barrier=(4.79431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87133,'amu*angstrom^2'), symmetry=1, barrier=(66.0176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4073.83,'J/mol'), sigma=(6.68315,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.32 K, Pc=30.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772443,0.0746308,-9.04976e-05,6.40734e-08,-1.87126e-11,82617.6,30.9806], Tmin=(100,'K'), Tmax=(829.3,'K')), NASAPolynomial(coeffs=[9.53949,0.0323443,-1.40115e-05,2.58676e-09,-1.76837e-13,81163.5,-9.6744], Tmin=(829.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C[C]=[C]CO[O](20095)',
    structure = SMILES('C=[C]C[C]=[C]CO[O]'),
    E0 = (836.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.237036,'amu*angstrom^2'), symmetry=1, barrier=(5.44993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236946,'amu*angstrom^2'), symmetry=1, barrier=(5.44785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237041,'amu*angstrom^2'), symmetry=1, barrier=(5.45003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236793,'amu*angstrom^2'), symmetry=1, barrier=(5.44434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59633,0.0845956,-0.000135462,1.23326e-07,-4.37326e-11,100736,33.5417], Tmin=(100,'K'), Tmax=(849.178,'K')), NASAPolynomial(coeffs=[6.15037,0.0383385,-1.82557e-05,3.44404e-09,-2.34792e-13,100517,11.9208], Tmin=(849.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(836.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C[C]=[C]COO(20324)',
    structure = SMILES('[CH]=[C]C[C]=[C]COO'),
    E0 = (931.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1685,1700,300,370,440,3120,650,792.5,1650,207.166,207.168],'cm^-1')),
        HinderedRotor(inertia=(0.316926,'amu*angstrom^2'), symmetry=1, barrier=(9.65211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316924,'amu*angstrom^2'), symmetry=1, barrier=(9.65209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00392789,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316928,'amu*angstrom^2'), symmetry=1, barrier=(9.65211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37269,'amu*angstrom^2'), symmetry=1, barrier=(41.8059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258273,0.0917151,-0.000147455,1.31545e-07,-4.60254e-11,112185,34.3994], Tmin=(100,'K'), Tmax=(834.469,'K')), NASAPolynomial(coeffs=[8.07826,0.0369957,-1.81138e-05,3.4625e-09,-2.38095e-13,111480,1.68301], Tmin=(834.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(931.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCC=[C][CH]O[O](20325)',
    structure = SMILES('[CH]=CCC=[C][CH]O[O]'),
    E0 = (725.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793782,0.0678201,-6.18318e-05,2.93149e-08,-5.62335e-12,87354.8,32.8886], Tmin=(100,'K'), Tmax=(1243.49,'K')), NASAPolynomial(coeffs=[13.8988,0.0256654,-1.09824e-05,2.05377e-09,-1.42694e-13,84095.6,-33.1916], Tmin=(1243.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(725.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]CC=[C][CH]O[O](20094)',
    structure = SMILES('C=[C]CC=[C][CH]O[O]'),
    E0 = (716.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,377.392,377.404,377.425],'cm^-1')),
        HinderedRotor(inertia=(0.00118376,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106928,'amu*angstrom^2'), symmetry=1, barrier=(10.8053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106925,'amu*angstrom^2'), symmetry=1, barrier=(10.8054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30792,'amu*angstrom^2'), symmetry=1, barrier=(31.1202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00777,0.0668246,-6.26284e-05,3.16609e-08,-6.6203e-12,86231,32.1941], Tmin=(100,'K'), Tmax=(1130.98,'K')), NASAPolynomial(coeffs=[11.5519,0.0295327,-1.31687e-05,2.50642e-09,-1.75776e-13,83846,-19.9728], Tmin=(1130.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(716.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH]C=[C]COO(20283)',
    structure = SMILES('[CH]=C=C[CH][C]COO'),
    E0 = (755.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.317484,0.0837264,-9.98763e-05,6.40855e-08,-1.65654e-11,91027.1,29.3008], Tmin=(100,'K'), Tmax=(940.066,'K')), NASAPolynomial(coeffs=[13.1243,0.0292344,-1.29292e-05,2.42673e-09,-1.68299e-13,88619.2,-31.6936], Tmin=(940.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
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
    E0 = (845.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1091.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (845.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1091.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1332.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1393.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1368.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (854.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (854.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (977.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (909.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (924.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (870.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (884.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1006.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1001.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (932.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (880.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (845.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (918.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1181.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1119.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1174.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1295.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (964.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (985.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1079.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1010.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (987.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (996.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (991.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1054.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1108.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (910.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (940.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (918.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['C3H3(5450)', 'C#CCO[O](16808)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=C([CH2])C=[C]CO[O](19498)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['O2(S)(666)', '[CH]=[C]CC=C=C(16903)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]CC=[C]C[O](20309)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=[CH](9646)', '[CH2]C=[C]CO[O](19545)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]O[O](46)', '[C]=CC[C]=[CH](17531)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[C]=[C]CC=[C]CO[O](20310)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=[C]CC=C1COO1(20311)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=C1CC=C1CO[O](20312)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriND_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['O(4)', '[CH]=[C]CC=C1CO1(20313)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO_S;Y_rad_intra;OOJ] for rate rule [R2OO_S;Cd_rad_out;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=C=CC=CCO[O](20314)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=[C]CC=C=COO(20315)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=CCC=C=CO[O](19507)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'C#C[CH]C=[C]CO[O](20280)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH]=[C]CC=C=CO[O](20316)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]=[C]CC#CCO[O](20317)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C3H3(5450)', '[CH]=[C]CO[O](16805)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.0130473,'m^3/(mol*s)'), n=2.32138, Ea=(27.0133,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C][CH2](16918)', 'C#CCO[O](16808)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.109931,'m^3/(mol*s)'), n=2.3439, Ea=(23.3771,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cs;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O2(2)', '[CH]=[C]CC=C=C(16903)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.973e+06,'cm^3/(mol*s)'), n=2.037, Ea=(140.233,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Ca;OJ] for rate rule [Cds-HH_Ca;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 140.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['O2(2)', '[CH]=[C]CC=[C][CH2](16907)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C][CH2](16918)', '[CH]=[C]CO[O](16805)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH][C]=CC=[C]CO[O](20256)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH]=[C]CC=[C][CH]O[O](20318)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH]=[C]C[C]=[C]CO[O](20319)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=C[CH]C=[C]CO[O](20281)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1661e+10,'s^-1'), n=0.959259, Ea=(118.299,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=[C]CC=C[CH]O[O](20320)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]C[C]=CCO[O](20321)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=[C]CC=[C][CH]OO(20322)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.66219e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=CC[C]=[C]CO[O](20323)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=[C][CH]C=CCO[O](20282)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.9395e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['C=[C]C[C]=[C]CO[O](20095)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=[C]C[C]=[C]COO(20324)'],
    products = ['[CH]=[C]CC=[C]CO[O](19500)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.44313e+09,'s^-1'), n=0.985167, Ea=(177.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=CCC=[C][CH]O[O](20325)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(26847.4,'s^-1'), n=2.09453, Ea=(64.7206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5HJ_3;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['C=[C]CC=[C][CH]O[O](20094)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.70099e+07,'s^-1'), n=1.485, Ea=(94.7468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R6Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH]=[C][CH]C=[C]COO(20283)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.48874e+06,'s^-1'), n=0.972854, Ea=(72.9565,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;O_rad_out;XH_out] for rate rule [R6HJ_3;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4521',
    isomers = [
        '[CH]=[C]CC=[C]CO[O](19500)',
    ],
    reactants = [
        ('C3H3(5450)', 'C#CCO[O](16808)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4521',
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


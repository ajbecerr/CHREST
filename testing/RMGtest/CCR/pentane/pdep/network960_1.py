species(
    label = '[C-]=[O+]CC=CC(3644)',
    structure = SMILES('[C-]=[O+]CC=CC'),
    E0 = (319.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,288,311.326,3986.38,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.823214,'amu*angstrom^2'), symmetry=1, barrier=(19.2917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823214,'amu*angstrom^2'), symmetry=1, barrier=(19.2917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823214,'amu*angstrom^2'), symmetry=1, barrier=(19.2917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02982,0.050076,-4.04504e-05,1.86181e-08,-3.95728e-12,38547.7,-0.70319], Tmin=(100,'K'), Tmax=(1021.49,'K')), NASAPolynomial(coeffs=[6.06468,0.0342761,-1.7249e-05,3.47589e-09,-2.51357e-13,37723.4,-20.2549], Tmin=(1021.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsCs) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CO(2039)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(815.652,'J/mol'), sigma=(3.65,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.95,'angstroms^3'), rotrelaxcollnum=1.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'butene1(35)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,247.029],'cm^-1')),
        HinderedRotor(inertia=(0.178818,'amu*angstrom^2'), symmetry=1, barrier=(7.72187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177126,'amu*angstrom^2'), symmetry=1, barrier=(7.72778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58774,0.0232777,1.93416e-05,-3.55502e-08,1.36908e-11,-1918.73,14.575], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20513,0.0236362,-9.03154e-06,1.65394e-09,-1.1602e-13,-3797.32,-12.4424], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[C-]=[O+]CC=C(3637)',
    structure = SMILES('[C-]=[O+]CC=C'),
    E0 = (355.987,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,262.591,262.964,264.86,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.193419,'amu*angstrom^2'), symmetry=1, barrier=(9.42786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169229,'amu*angstrom^2'), symmetry=1, barrier=(9.37212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2836.07,'J/mol'), sigma=(5.784e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6126,0.0356299,-2.5255e-05,8.59323e-09,-1.22402e-12,42860.3,-4.63114], Tmin=(100,'K'), Tmax=(1507.78,'K')), NASAPolynomial(coeffs=[8.15635,0.0209232,-1.06244e-05,2.12444e-09,-1.51476e-13,41188.5,-33.6532], Tmin=(1507.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[C-][O+]=C(3687)',
    structure = SMILES('[C-][O+]=C'),
    E0 = (339.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36274,0.0208744,-4.99693e-05,5.9339e-08,-2.45823e-11,40899.7,-11.1374], Tmin=(100,'K'), Tmax=(842.434,'K')), NASAPolynomial(coeffs=[0.855816,0.0151766,-8.4845e-06,1.70879e-09,-1.20102e-13,41946.7,4.23414], Tmin=(842.434,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsJ2_singlet-CsH) + group(Cds-CdsCsCs)"""),
)

species(
    label = 'propyne(3682)',
    structure = SMILES('C#CC'),
    E0 = (172.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30523,0.0109265,1.31989e-05,-2.25163e-08,8.87417e-12,20738.2,7.19163], Tmin=(100,'K'), Tmax=(969.986,'K')), NASAPolynomial(coeffs=[5.80055,0.0113536,-4.03491e-06,7.19063e-10,-5.02185e-14,19750,-7.36966], Tmin=(969.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""propyne""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[C-]=[O+]C[C]CC(3692)',
    structure = SMILES('[C-]=[O+]C[C]CC'),
    E0 = (634.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,709.046,1600,2155.62,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0904173,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904173,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904173,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904173,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95855,0.055193,-2.19507e-05,-7.53831e-08,9.06399e-11,76342.2,9.65582], Tmin=(100,'K'), Tmax=(446.978,'K')), NASAPolynomial(coeffs=[5.77174,0.0377185,-1.9183e-05,3.82524e-09,-2.73173e-13,75835,-7.53061], Tmin=(446.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(634.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsCsCsCs) + group(CsJ2_singlet-CsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[C-]=[O+]CC[C]C(3693)',
    structure = SMILES('[C-]=[O+]CC[C]C'),
    E0 = (634.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,709.048,1600,2155.61,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0904172,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904172,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904172,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0904172,'amu*angstrom^2'), symmetry=1, barrier=(2.07887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95756,0.0552105,-2.20557e-05,-7.513e-08,9.04283e-11,76342.2,9.65908], Tmin=(100,'K'), Tmax=(447.056,'K')), NASAPolynomial(coeffs=[5.77189,0.0377182,-1.91828e-05,3.82518e-09,-2.73168e-13,75834.9,-7.5314], Tmin=(447.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(634.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsCsCsCs) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (395.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (774.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (616.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (671.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (671.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[C-]=[O+]CC=CC(3644)'],
    products = ['CO(2039)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.60496e+13,'s^-1'), n=-0.343905, Ea=(76.0365,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=6.4548336266697035, var=111.97246879746811, Tref=1000.0, N=4, correlation='Root_1R!H->C_N-2R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_N-2R!H->C
    Total Standard Deviation in ln(k): 37.43168900536777
Exact match found for rate rule [Root_1R!H->C_N-2R!H->C]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[C-]=[O+]CC=C(3637)'],
    products = ['[C-]=[O+]CC=CC(3644)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C-]=[O+]CC=CC(3644)'],
    products = ['[C-][O+]=C(3687)', 'propyne(3682)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.29914e+17,'s^-1'), n=-1.73308, Ea=(296.513,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.16692063000437535, var=8.814681701569807, Tref=1000.0, N=65, correlation='Root',), comment="""BM rule fitted to 2 training reactions at node Root
    Total Standard Deviation in ln(k): 6.371362717619001
Exact match found for rate rule [Root]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C-]=[O+]C[C]CC(3692)'],
    products = ['[C-]=[O+]CC=CC(3644)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C-]=[O+]CC[C]C(3693)'],
    products = ['[C-]=[O+]CC=CC(3644)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = 'PDepNetwork #960',
    isomers = [
        '[C-]=[O+]CC=CC(3644)',
    ],
    reactants = [
        ('CO(2039)', 'butene1(35)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #960',
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


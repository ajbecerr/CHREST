species(
    label = '[CH2]CO[C]=C[O](10945)',
    structure = SMILES('[CH2]CO[C]=C[O]'),
    E0 = (198.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,264.998,265.638,265.684,265.957],'cm^-1')),
        HinderedRotor(inertia=(0.427169,'amu*angstrom^2'), symmetry=1, barrier=(21.5405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428673,'amu*angstrom^2'), symmetry=1, barrier=(21.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.430742,'amu*angstrom^2'), symmetry=1, barrier=(21.5346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985411,0.0533117,-2.46087e-05,-2.27174e-08,1.74139e-11,23994.2,24.5369], Tmin=(100,'K'), Tmax=(914.987,'K')), NASAPolynomial(coeffs=[20.4274,0.00288864,1.37984e-06,-3.60174e-10,2.29684e-14,18989.3,-75.4398], Tmin=(914.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H4(22)',
    structure = SMILES('C=C'),
    E0 = (41.9072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98876,-0.00674712,5.04398e-05,-5.70744e-08,2.04945e-11,5047.05,3.80494], Tmin=(100,'K'), Tmax=(946.013,'K')), NASAPolynomial(coeffs=[4.59022,0.00872722,-2.66493e-06,4.8171e-10,-3.6069e-14,4127.02,-3.3246], Tmin=(946.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.9072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[C]=C[O](6859)',
    structure = SMILES('[C]C=O'),
    E0 = (491.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53714,0.00730776,2.32451e-06,-8.52039e-09,3.75032e-12,59149.2,7.56143], Tmin=(100,'K'), Tmax=(1018.09,'K')), NASAPolynomial(coeffs=[6.47558,0.00262438,-8.8465e-07,2.00876e-10,-1.68063e-14,58195.3,-8.41394], Tmin=(1018.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJ3)"""),
)

species(
    label = '[CH2]C[O](1358)',
    structure = SMILES('[CH2]C[O]'),
    E0 = (187.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1795.79],'cm^-1')),
        HinderedRotor(inertia=(0.00405474,'amu*angstrom^2'), symmetry=1, barrier=(9.27899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3197.09,'J/mol'), sigma=(5.50868,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=499.38 K, Pc=43.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63182,0.0106853,2.91824e-06,-4.98315e-09,1.21781e-12,22581.3,11.2346], Tmin=(100,'K'), Tmax=(1659.64,'K')), NASAPolynomial(coeffs=[5.14351,0.0138609,-6.11494e-06,1.12108e-09,-7.46163e-14,21140.5,0.346415], Tmin=(1659.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO)"""),
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
    label = '[CH2]O[C]=C[O](11426)',
    structure = SMILES('[CH2]O[C]=C[O]'),
    E0 = (216.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,400.773,400.953,401.443],'cm^-1')),
        HinderedRotor(inertia=(0.189956,'amu*angstrom^2'), symmetry=1, barrier=(21.6878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18896,'amu*angstrom^2'), symmetry=1, barrier=(21.685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09104,0.028776,1.00821e-05,-4.80262e-08,2.47341e-11,26115.5,20.59], Tmin=(100,'K'), Tmax=(909.214,'K')), NASAPolynomial(coeffs=[17.1662,-0.00241134,3.5706e-06,-7.50851e-10,4.89516e-14,21922,-58.6903], Tmin=(909.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COCJ) + radical(C=CJO)"""),
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
    label = '[CH]=[C]OC[CH2](2368)',
    structure = SMILES('[CH]=[C]OC[CH2]'),
    E0 = (512.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,351.168,351.747],'cm^-1')),
        HinderedRotor(inertia=(0.191788,'amu*angstrom^2'), symmetry=1, barrier=(16.7995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191346,'amu*angstrom^2'), symmetry=1, barrier=(16.8045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191479,'amu*angstrom^2'), symmetry=1, barrier=(16.806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3055.61,'J/mol'), sigma=(5.41904,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=477.28 K, Pc=43.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19429,0.0520531,-5.40358e-05,2.74157e-08,-5.27588e-12,61799.6,22.9226], Tmin=(100,'K'), Tmax=(1404.71,'K')), NASAPolynomial(coeffs=[15.1524,0.0076461,-1.6399e-06,1.8709e-10,-9.59633e-15,58337.9,-47.5241], Tmin=(1404.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CJCO) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]CO[C]=C[O](11427)',
    structure = SMILES('[CH]CO[C]=C[O]'),
    E0 = (435.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08917,0.0511297,-2.41925e-05,-2.24588e-08,1.7316e-11,52450,24.5349], Tmin=(100,'K'), Tmax=(913.759,'K')), NASAPolynomial(coeffs=[20.487,0.000630124,2.21168e-06,-5.05313e-10,3.26694e-14,47468.2,-75.1606], Tmin=(913.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]COC1=CO1(11428)',
    structure = SMILES('[CH2]COC1=CO1'),
    E0 = (191.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533931,0.0639664,-6.80037e-05,3.45416e-08,-6.61836e-12,23143.8,21.3014], Tmin=(100,'K'), Tmax=(1410.81,'K')), NASAPolynomial(coeffs=[18.5584,0.00707343,-1.35879e-06,1.40539e-10,-6.96159e-15,18634.1,-69.8178], Tmin=(1410.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(CJCO)"""),
)

species(
    label = '[O]C=C1CCO1(11429)',
    structure = SMILES('[O]C=C1CCO1'),
    E0 = (-122.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83359,0.0223984,7.47765e-05,-1.35441e-07,6.09198e-11,-14639.6,15.3507], Tmin=(100,'K'), Tmax=(895.445,'K')), NASAPolynomial(coeffs=[23.1436,-0.00534429,8.2605e-06,-1.79816e-09,1.22033e-13,-21160.2,-100.204], Tmin=(895.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ)"""),
)

species(
    label = '[CH2]COC=C=O(2796)',
    structure = SMILES('[CH2]COC=C=O'),
    E0 = (-1.23474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341696,0.0636452,-6.66348e-05,3.29358e-08,-6.0394e-12,-2.36848,23.121], Tmin=(100,'K'), Tmax=(1544.5,'K')), NASAPolynomial(coeffs=[18.9333,0.00522996,6.76534e-08,-1.59134e-10,1.40783e-14,-4520.84,-70.6906], Tmin=(1544.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1.23474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CJCO)"""),
)

species(
    label = 'C=COC=C[O](2765)',
    structure = SMILES('C=COC=C[O]'),
    E0 = (-97.5162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38427,0.0384861,2.24397e-05,-7.41186e-08,3.63062e-11,-11616.6,21.6985], Tmin=(100,'K'), Tmax=(916.099,'K')), NASAPolynomial(coeffs=[21.5709,0.000258613,3.304e-06,-7.17302e-10,4.45496e-14,-17409.7,-85.3531], Tmin=(916.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.5162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C=CO[C]=C[O](11102)',
    structure = SMILES('C=CO[C]=C[O]'),
    E0 = (142.228,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12889,'amu*angstrom^2'), symmetry=1, barrier=(25.9555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12474,'amu*angstrom^2'), symmetry=1, barrier=(25.8601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59936,0.0401881,-4.28583e-06,-3.54569e-08,2.00889e-11,17204.2,24.0548], Tmin=(100,'K'), Tmax=(921.403,'K')), NASAPolynomial(coeffs=[17.409,0.00442144,4.35709e-07,-1.6052e-10,8.40455e-15,12895.7,-58.4943], Tmin=(921.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CO[C]=C=O(11430)',
    structure = SMILES('[CH2]CO[C]=C=O'),
    E0 = (238.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,381.232,381.234,381.252],'cm^-1')),
        HinderedRotor(inertia=(0.146581,'amu*angstrom^2'), symmetry=1, barrier=(15.1209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146623,'amu*angstrom^2'), symmetry=1, barrier=(15.1214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146593,'amu*angstrom^2'), symmetry=1, barrier=(15.1212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18828,0.0580469,-6.84278e-05,3.99546e-08,-8.99337e-12,28790.6,23.201], Tmin=(100,'K'), Tmax=(1099.38,'K')), NASAPolynomial(coeffs=[14.2719,0.0104431,-3.47674e-06,5.68057e-10,-3.6804e-14,25913.8,-41.1595], Tmin=(1099.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = 'C#COC[CH2](2366)',
    structure = SMILES('C#COC[CH2]'),
    E0 = (235.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17084,'amu*angstrom^2'), symmetry=1, barrier=(26.9199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17373,'amu*angstrom^2'), symmetry=1, barrier=(26.9864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17381,'amu*angstrom^2'), symmetry=1, barrier=(26.9881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789011,0.0554533,-5.66782e-05,2.73257e-08,-4.90149e-12,28441.7,18.3114], Tmin=(100,'K'), Tmax=(1571.19,'K')), NASAPolynomial(coeffs=[17.2644,0.00496135,-3.13134e-07,-5.30374e-11,5.8481e-15,24319.7,-65.2593], Tmin=(1571.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCO)"""),
)

species(
    label = 'C2H4(T)(23)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,995.342,995.492,2246.13,2948.67],'cm^-1')),
        HinderedRotor(inertia=(0.00723049,'amu*angstrom^2'), symmetry=1, barrier=(5.08405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40735,0.0100313,6.40897e-06,-1.41287e-08,5.92653e-12,38288.2,6.11706], Tmin=(100,'K'), Tmax=(954.266,'K')), NASAPolynomial(coeffs=[5.52251,0.00856169,-2.90741e-06,5.02347e-10,-3.44567e-14,37547.7,-5.75289], Tmin=(954.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]O[C]=C[O](10214)',
    structure = SMILES('[CH2][CH]O[C]=C[O]'),
    E0 = (392.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.967915,'amu*angstrom^2'), symmetry=1, barrier=(22.2543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.969452,'amu*angstrom^2'), symmetry=1, barrier=(22.2896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968215,'amu*angstrom^2'), symmetry=1, barrier=(22.2612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.246703,0.0728801,-8.70227e-05,4.61145e-08,-8.82495e-12,47368.2,29.1271], Tmin=(100,'K'), Tmax=(1521.7,'K')), NASAPolynomial(coeffs=[22.602,-0.00350478,4.3639e-06,-9.72175e-10,6.90951e-14,42304.4,-84.4874], Tmin=(1521.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CO[C]=[C][O](10215)',
    structure = SMILES('[CH2]CO[C]=[C][O]'),
    E0 = (438.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,432.173,432.174,432.175,432.178],'cm^-1')),
        HinderedRotor(inertia=(0.126195,'amu*angstrom^2'), symmetry=1, barrier=(16.7259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126196,'amu*angstrom^2'), symmetry=1, barrier=(16.7259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126196,'amu*angstrom^2'), symmetry=1, barrier=(16.7259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.838684,0.0592763,-6.62264e-05,3.52495e-08,-7.04398e-12,52830.8,28.1914], Tmin=(100,'K'), Tmax=(1363.81,'K')), NASAPolynomial(coeffs=[17.129,0.00556422,-6.25031e-07,-8.10147e-12,3.8387e-15,48939.2,-53.4315], Tmin=(1363.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]O[C]=C[O](10373)',
    structure = SMILES('C[CH]O[C]=C[O]'),
    E0 = (180.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.959993,'amu*angstrom^2'), symmetry=1, barrier=(22.0721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958195,'amu*angstrom^2'), symmetry=1, barrier=(22.0308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958025,'amu*angstrom^2'), symmetry=1, barrier=(22.0269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.372209,0.0707098,-7.67305e-05,3.78118e-08,-6.7487e-12,21928.8,28.2688], Tmin=(100,'K'), Tmax=(1648.76,'K')), NASAPolynomial(coeffs=[21.4681,0.000503878,2.80736e-06,-6.83452e-10,4.88076e-14,17067.4,-80.9208], Tmin=(1648.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CO[C]=[C]O(11431)',
    structure = SMILES('[CH2]CO[C]=[C]O'),
    E0 = (296.778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,286.157,286.185,286.239],'cm^-1')),
        HinderedRotor(inertia=(0.293154,'amu*angstrom^2'), symmetry=1, barrier=(17.0565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293542,'amu*angstrom^2'), symmetry=1, barrier=(17.0575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29332,'amu*angstrom^2'), symmetry=1, barrier=(17.056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293702,'amu*angstrom^2'), symmetry=1, barrier=(17.0568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.178228,0.0692284,-8.0307e-05,4.32943e-08,-8.58376e-12,35844.5,29.1487], Tmin=(100,'K'), Tmax=(1435.1,'K')), NASAPolynomial(coeffs=[19.9062,0.00250148,1.70821e-06,-5.05605e-10,3.93568e-14,31391.1,-68.942], Tmin=(1435.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CO[CH][C]=O(2798)',
    structure = SMILES('[CH2]COC=[C][O]'),
    E0 = (198.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,264.998,265.638,265.684,265.957],'cm^-1')),
        HinderedRotor(inertia=(0.427169,'amu*angstrom^2'), symmetry=1, barrier=(21.5405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428673,'amu*angstrom^2'), symmetry=1, barrier=(21.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.430742,'amu*angstrom^2'), symmetry=1, barrier=(21.5346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985411,0.0533117,-2.46087e-05,-2.27174e-08,1.74139e-11,23994.2,24.5369], Tmin=(100,'K'), Tmax=(914.987,'K')), NASAPolynomial(coeffs=[20.4274,0.00288864,1.37984e-06,-3.60174e-10,2.29684e-14,18989.3,-75.4398], Tmin=(914.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]OC=C[O](1927)',
    structure = SMILES('[CH2][CH]OC=C[O]'),
    E0 = (152.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13522,'amu*angstrom^2'), symmetry=1, barrier=(26.101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13573,'amu*angstrom^2'), symmetry=1, barrier=(26.1127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13566,'amu*angstrom^2'), symmetry=1, barrier=(26.1111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28393,0.080444,-9.0962e-05,4.51811e-08,-7.97971e-12,18584.4,29.7525], Tmin=(100,'K'), Tmax=(1707.69,'K')), NASAPolynomial(coeffs=[24.2187,-0.00453177,5.84919e-06,-1.26816e-09,8.7385e-14,13554.5,-96.155], Tmin=(1707.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CJCO)"""),
)

species(
    label = 'CCO[C]=[C][O](10374)',
    structure = SMILES('CCO[C]=[C][O]'),
    E0 = (226.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,390.522,390.6,390.762,391.278],'cm^-1')),
        HinderedRotor(inertia=(0.149455,'amu*angstrom^2'), symmetry=1, barrier=(16.2084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(16.2078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149455,'amu*angstrom^2'), symmetry=1, barrier=(16.2078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34198,0.0499585,-3.20844e-05,-2.52079e-09,7.03461e-12,27363.3,25.7511], Tmin=(100,'K'), Tmax=(941.543,'K')), NASAPolynomial(coeffs=[15.4991,0.0107505,-2.97545e-06,4.85329e-10,-3.43828e-14,23769.4,-46.6244], Tmin=(941.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]O[C]=CO(11432)',
    structure = SMILES('[CH2][CH]O[C]=CO'),
    E0 = (250.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.936266,'amu*angstrom^2'), symmetry=1, barrier=(21.5266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.935979,'amu*angstrom^2'), symmetry=1, barrier=(21.52,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934266,'amu*angstrom^2'), symmetry=1, barrier=(21.4806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940648,'amu*angstrom^2'), symmetry=1, barrier=(21.6273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.912758,0.0828938,-0.000101299,5.43865e-08,-1.04498e-11,30382.2,30.1049], Tmin=(100,'K'), Tmax=(1537.76,'K')), NASAPolynomial(coeffs=[25.0146,-0.00605176,6.43814e-06,-1.41461e-09,1.00402e-13,24950.7,-97.87], Tmin=(1537.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOC(O)) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[C]1[CH]OCCO1(11433)',
    structure = SMILES('[C]1[CH]OCCO1'),
    E0 = (302.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355918,0.0583938,-6.03458e-05,3.32704e-08,-6.54146e-12,36527.6,18.6255], Tmin=(100,'K'), Tmax=(1610.55,'K')), NASAPolynomial(coeffs=[9.5061,0.0136897,1.75949e-06,-9.1043e-10,7.95713e-14,36430.7,-21.0302], Tmin=(1610.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(1,4-Dioxane) + radical(CCsJOCs) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1[C]OCC1(11434)',
    structure = SMILES('[O]C1[C]OCC1'),
    E0 = (320.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19049,0.0261537,3.66218e-05,-7.57657e-08,3.49706e-11,38603.3,14.9139], Tmin=(100,'K'), Tmax=(882.858,'K')), NASAPolynomial(coeffs=[13.9577,0.0102148,2.00705e-07,-3.11563e-10,2.56373e-14,35069,-48.6396], Tmin=(882.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CH2_triplet)"""),
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
    E0 = (198.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (679.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (654.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1032.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (646.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (201.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (206.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (220.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (287.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (360.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (458.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (500.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (272.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (522.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (604.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (650.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (356.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (455.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (368.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (340.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (259.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (303.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (302.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (320.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['OCHCO(3676)', 'C2H4(22)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=C[O](6859)', '[CH2]C[O](1358)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(20)', '[CH2]O[C]=C[O](11426)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]OC[CH2](2368)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]CO[C]=C[O](11427)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['[CH2]COC1=CO1(11428)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['[O]C=C1CCO1(11429)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['[CH2]COC=C=O(2796)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['C=COC=C[O](2765)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C=CO[C]=C[O](11102)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]CO[C]=C=O(11430)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O(4)', 'C#COC[CH2](2366)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.30864e+06,'m^3/(mol*s)'), n=-0.19959, Ea=(22.3126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_atom_triplet] + [Ct_Ct;YJ] for rate rule [Ct_Ct;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OCHCO(3676)', 'C2H4(T)(23)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(23.3993,'m^3/(mol*s)'), n=2.021, Ea=(29.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=C[O](9592)', 'C2H4(T)(23)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.62e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][CH]O[C]=C[O](10214)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]CO[C]=[C][O](10215)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['C[CH]O[C]=C[O](10373)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.7e+13,'s^-1','+|-',2), n=-0.1, Ea=(158.364,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""From training reaction 347 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CO[C]=[C]O(11431)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CO[CH][C]=O(2798)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['[CH2][CH]OC=C[O](1927)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCO[C]=[C][O](10374)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]O[C]=CO(11432)'],
    products = ['[CH2]CO[C]=C[O](10945)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.10205e+06,'s^-1'), n=1.54368, Ea=(52.1315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;XH_out] for rate rule [R5HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['[C]1[CH]OCCO1(11433)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(103.964,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 102.7 to 104.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CO[C]=C[O](10945)'],
    products = ['[O]C1[C]OCC1(11434)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.93455e+07,'s^-1'), n=0.812267, Ea=(121.82,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 118.1 to 121.8 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #2591',
    isomers = [
        '[CH2]CO[C]=C[O](10945)',
    ],
    reactants = [
        ('OCHCO(3676)', 'C2H4(22)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2591',
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


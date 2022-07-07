species(
    label = '[CH]=C(C=[C][O])CO[O](19438)',
    structure = SMILES('[CH]C(=C[C]=O)CO[O]'),
    E0 = (445.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.21403,'amu*angstrom^2'), symmetry=1, barrier=(50.9048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20759,'amu*angstrom^2'), symmetry=1, barrier=(50.7568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21129,'amu*angstrom^2'), symmetry=1, barrier=(50.8418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20475,'amu*angstrom^2'), symmetry=1, barrier=(50.6914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.226244,0.0939598,-0.000156483,1.43949e-07,-5.16224e-11,53677.8,30.4513], Tmin=(100,'K'), Tmax=(832.305,'K')), NASAPolynomial(coeffs=[7.03552,0.0398085,-2.02754e-05,3.91814e-09,-2.70667e-13,53286.4,3.30871], Tmin=(832.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = '[CH]C(=C[C]=O)C[O](20739)',
    structure = SMILES('[CH]C(=C[C]=O)C[O]'),
    E0 = (447.456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,361.99,363.236,363.307,363.982,364.151],'cm^-1')),
        HinderedRotor(inertia=(0.556509,'amu*angstrom^2'), symmetry=1, barrier=(52.016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551834,'amu*angstrom^2'), symmetry=1, barrier=(52.0018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555474,'amu*angstrom^2'), symmetry=1, barrier=(52.0057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998229,0.0779242,-0.000133215,1.31966e-07,-5.03669e-11,53913,27.0187], Tmin=(100,'K'), Tmax=(829.288,'K')), NASAPolynomial(coeffs=[2.74473,0.0430534,-2.23055e-05,4.34959e-09,-3.02218e-13,54532.8,24.4028], Tmin=(829.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[C]=O(2355)',
    structure = SMILES('[C]=O'),
    E0 = (438.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3323.79],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09152,0.00193295,-1.59153e-05,2.47563e-08,-1.11287e-11,52770.5,4.46624], Tmin=(100,'K'), Tmax=(866.029,'K')), NASAPolynomial(coeffs=[1.05092,0.00526657,-3.13864e-06,6.40624e-10,-4.47939e-14,53698.8,21.0169], Tmin=(866.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])CO[O](20211)',
    structure = SMILES('[CH]C(=[CH])CO[O]'),
    E0 = (655.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,365.959,366.184,366.425],'cm^-1')),
        HinderedRotor(inertia=(0.55896,'amu*angstrom^2'), symmetry=1, barrier=(53.0299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.55651,'amu*angstrom^2'), symmetry=1, barrier=(53.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555869,'amu*angstrom^2'), symmetry=1, barrier=(53.0372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52129,0.0585443,-7.59095e-05,6.10793e-08,-2.03843e-11,78890.5,22.9153], Tmin=(100,'K'), Tmax=(819.602,'K')), NASAPolynomial(coeffs=[6.42814,0.0294529,-1.32536e-05,2.45733e-09,-1.67299e-13,78259,1.27275], Tmin=(819.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1=CC(=O)OOC1(20740)',
    structure = SMILES('[CH]=C1[CH]C(=O)OOC1'),
    E0 = (59.3471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07008,0.0795541,-5.89274e-05,1.51897e-08,-1.4074e-12,7169.17,0.0488668], Tmin=(100,'K'), Tmax=(2597.56,'K')), NASAPolynomial(coeffs=[58.0499,-0.00665131,-9.14572e-06,2.41296e-09,-1.77687e-13,-21912.4,-323.457], Tmin=(2597.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.3471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(six-sidedoubles) + radical(C=CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1=CC(=O)OC1(20741)',
    structure = SMILES('[CH]=C1[CH]C(=O)OC1'),
    E0 = (76.9004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.583,0.00913788,9.06715e-05,-1.21719e-07,4.54116e-11,9319.61,17.9052], Tmin=(100,'K'), Tmax=(990.796,'K')), NASAPolynomial(coeffs=[14.2024,0.0182594,-7.9646e-06,1.72528e-09,-1.37903e-13,4266.92,-51.9228], Tmin=(990.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.9004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C=C=O)COO(20742)',
    structure = SMILES('[CH]=C([C]=C=O)COO'),
    E0 = (339.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646475,0.081771,-0.000131737,1.16734e-07,-4.07491e-11,40952.2,29.9195], Tmin=(100,'K'), Tmax=(822.702,'K')), NASAPolynomial(coeffs=[8.24105,0.0315039,-1.57615e-05,3.04242e-09,-2.10549e-13,40154.1,-2.49364], Tmin=(822.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)CJ=C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]C([CH][C]=O)[CH]O[O](20743)',
    structure = SMILES('[CH]C([CH][C]=O)[CH]O[O]'),
    E0 = (768.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366225,0.0872872,-0.000141357,1.19657e-07,-3.90543e-11,92498.3,32.6005], Tmin=(100,'K'), Tmax=(879.48,'K')), NASAPolynomial(coeffs=[10.4479,0.0264281,-1.19644e-05,2.17408e-09,-1.43941e-13,91305.3,-11.4438], Tmin=(879.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(768.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C1(CO[O])[CH]C1=O(20744)',
    structure = SMILES('[CH]C1(CO[O])[CH]C1=O'),
    E0 = (610.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.206052,0.0737884,-7.91664e-05,4.17655e-08,-8.47067e-12,73550.8,28.3488], Tmin=(100,'K'), Tmax=(1248.3,'K')), NASAPolynomial(coeffs=[18.8669,0.0125367,-3.81503e-06,5.8919e-10,-3.70872e-14,69005.4,-65.3632], Tmin=(1248.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(ROOJ) + radical(CCJC=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1COOC1[C]=O(20745)',
    structure = SMILES('[CH][C]1COOC1[C]=O'),
    E0 = (525.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8126,0.0287213,5.28727e-05,-1.0327e-07,4.6811e-11,63265.3,31.132], Tmin=(100,'K'), Tmax=(895.023,'K')), NASAPolynomial(coeffs=[18.2726,0.00707814,2.13206e-06,-6.62679e-10,4.68687e-14,58239.3,-58.0697], Tmin=(895.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(C2CJCOOH) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH][C]=O)COO1(20746)',
    structure = SMILES('[CH]C1([CH][C]=O)COO1'),
    E0 = (607.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.765261,0.075468,-0.000103934,7.82184e-08,-2.36785e-11,73189.9,27.2044], Tmin=(100,'K'), Tmax=(807.376,'K')), NASAPolynomial(coeffs=[10.7579,0.0259581,-1.19452e-05,2.25645e-09,-1.55778e-13,71576.4,-18.8656], Tmin=(807.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(CCJCO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = '[CH]C(=C=C=O)CO[O](20747)',
    structure = SMILES('[CH]=C([C]=C=O)CO[O]'),
    E0 = (491.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,3120,650,792.5,1650,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.315727,'amu*angstrom^2'), symmetry=1, barrier=(7.2592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315373,'amu*angstrom^2'), symmetry=1, barrier=(7.25104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315899,'amu*angstrom^2'), symmetry=1, barrier=(7.26313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.88941,0.0784232,-0.000137233,1.26853e-07,-4.46884e-11,59223.5,29.8063], Tmin=(100,'K'), Tmax=(861.173,'K')), NASAPolynomial(coeffs=[6.70088,0.0297584,-1.472e-05,2.78957e-09,-1.89401e-13,59026.2,7.30365], Tmin=(861.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=C(C)CJ=C=O) + radical(Cds_P)"""),
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
    label = '[CH]C([CH2])=C[C]=O(16420)',
    structure = SMILES('[CH]C([CH2])=C[C]=O'),
    E0 = (526.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10045,'amu*angstrom^2'), symmetry=1, barrier=(48.2934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09767,'amu*angstrom^2'), symmetry=1, barrier=(48.2295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10075,'amu*angstrom^2'), symmetry=1, barrier=(48.3004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.19,'J/mol'), sigma=(5.87487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.81 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52118,0.0593111,-6.47653e-05,4.05501e-08,-1.0837e-11,63394.7,20.8288], Tmin=(100,'K'), Tmax=(884.972,'K')), NASAPolynomial(coeffs=[8.03996,0.0298466,-1.48236e-05,2.92786e-09,-2.08896e-13,62241,-9.82393], Tmin=(884.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = '[CH][C]=C[C]=O(10445)',
    structure = SMILES('[CH][C]=C[C]=O'),
    E0 = (651.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07768,'amu*angstrom^2'), symmetry=1, barrier=(47.7699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08207,'amu*angstrom^2'), symmetry=1, barrier=(47.8709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.0581,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99066,0.0507484,-8.43786e-05,8.05615e-08,-3.01685e-11,78457.4,18.5596], Tmin=(100,'K'), Tmax=(815.064,'K')), NASAPolynomial(coeffs=[4.49774,0.025284,-1.3295e-05,2.60942e-09,-1.82246e-13,78485.9,9.65873], Tmin=(815.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH]O[O])=C[C]=O(20748)',
    structure = SMILES('[CH]C([CH]O[O])=C[C]=O'),
    E0 = (562.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,189.826,189.826,189.827,189.827],'cm^-1')),
        HinderedRotor(inertia=(1.93466,'amu*angstrom^2'), symmetry=1, barrier=(49.4702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93465,'amu*angstrom^2'), symmetry=1, barrier=(49.4702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93465,'amu*angstrom^2'), symmetry=1, barrier=(49.4702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93465,'amu*angstrom^2'), symmetry=1, barrier=(49.4702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639997,0.0801488,-0.0001052,7.60532e-08,-2.26019e-11,67775.2,29.5845], Tmin=(100,'K'), Tmax=(813.733,'K')), NASAPolynomial(coeffs=[10.4478,0.0319361,-1.63249e-05,3.23824e-09,-2.30665e-13,66179.1,-15.7106], Tmin=(813.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C(=[C][C]=O)CO[O](20749)',
    structure = SMILES('[CH]C(=[C][C]=O)CO[O]'),
    E0 = (689.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.12937,'amu*angstrom^2'), symmetry=1, barrier=(48.9584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13417,'amu*angstrom^2'), symmetry=1, barrier=(49.0687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13347,'amu*angstrom^2'), symmetry=1, barrier=(49.0527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1276,'amu*angstrom^2'), symmetry=1, barrier=(48.9176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0685086,0.0990788,-0.000177097,1.65136e-07,-5.8981e-11,83023.6,31.722], Tmin=(100,'K'), Tmax=(842.111,'K')), NASAPolynomial(coeffs=[7.85653,0.0362659,-1.92207e-05,3.74189e-09,-2.58297e-13,82627.5,0.923464], Tmin=(842.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C(=[C]C=O)CO[O](20750)',
    structure = SMILES('[CH]C(=C=C[O])CO[O]'),
    E0 = (481.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15859,'amu*angstrom^2'), symmetry=1, barrier=(49.6303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15796,'amu*angstrom^2'), symmetry=1, barrier=(49.6158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15779,'amu*angstrom^2'), symmetry=1, barrier=(49.6118,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431362,0.0781722,-8.66879e-05,5.08128e-08,-1.19255e-11,58025,29.249], Tmin=(100,'K'), Tmax=(1035.46,'K')), NASAPolynomial(coeffs=[13.984,0.0258175,-1.08443e-05,1.98139e-09,-1.35572e-13,55218.4,-36.6067], Tmin=(1035.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH]OO)=C[C]=O(20751)',
    structure = SMILES('[CH]C([CH]OO)=C[C]=O'),
    E0 = (410.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564797,0.0814756,-9.23577e-05,5.57691e-08,-1.3915e-11,49496.6,29.0983], Tmin=(100,'K'), Tmax=(956.735,'K')), NASAPolynomial(coeffs=[12.0133,0.0336102,-1.73122e-05,3.4759e-09,-2.50394e-13,47306,-25.6277], Tmin=(956.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJO) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C([CH]O[O])=CC=O(20752)',
    structure = SMILES('[CH]C([CH]O[O])=CC=O'),
    E0 = (401.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22964,0.0658954,-5.87376e-05,2.80314e-08,-5.72832e-12,48437.1,27.934], Tmin=(100,'K'), Tmax=(1121.3,'K')), NASAPolynomial(coeffs=[9.89556,0.034981,-1.73818e-05,3.443e-09,-2.46112e-13,46493.7,-14.866], Tmin=(1121.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[C][C]=O)COO(20753)',
    structure = SMILES('[CH]C(=[C][C]=O)COO'),
    E0 = (537.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.170599,0.102388,-0.000171516,1.55027e-07,-5.5129e-11,64752.1,31.8212], Tmin=(100,'K'), Tmax=(808.214,'K')), NASAPolynomial(coeffs=[9.34924,0.038095,-2.03116e-05,4.00664e-09,-2.80445e-13,63774.4,-8.60872], Tmin=(808.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=[C]C([C]=O)CO[O](20754)',
    structure = SMILES('[CH]=[C]C([C]=O)CO[O]'),
    E0 = (577.312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.351284,'amu*angstrom^2'), symmetry=1, barrier=(8.07671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.351118,'amu*angstrom^2'), symmetry=1, barrier=(8.0729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.351166,'amu*angstrom^2'), symmetry=1, barrier=(8.07399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.351219,'amu*angstrom^2'), symmetry=1, barrier=(8.07521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493896,0.0883413,-0.000156581,1.43951e-07,-5.01888e-11,69550.1,32.3897], Tmin=(100,'K'), Tmax=(868.835,'K')), NASAPolynomial(coeffs=[7.54815,0.0311945,-1.53285e-05,2.8854e-09,-1.94646e-13,69255.5,4.7074], Tmin=(868.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(CC(C)CJ=O) + radical(Cds_P)"""),
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
    label = '[C]=C([CH][C]=O)CO[O](20755)',
    structure = SMILES('[C]C(=C[C]=O)CO[O]'),
    E0 = (744.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.683439,'amu*angstrom^2'), symmetry=1, barrier=(15.7136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682352,'amu*angstrom^2'), symmetry=1, barrier=(15.6886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68258,'amu*angstrom^2'), symmetry=1, barrier=(15.6939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.150948,0.096413,-0.000177844,1.63767e-07,-5.73752e-11,89616.3,28.0659], Tmin=(100,'K'), Tmax=(846.093,'K')), NASAPolynomial(coeffs=[9.63126,0.0281725,-1.53406e-05,3.00807e-09,-2.07751e-13,88850.4,-11.1326], Tmin=(846.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJ=O) + radical(CJ3)"""),
)

species(
    label = '[O]OCC1=CC1[C]=O(20756)',
    structure = SMILES('[O]OCC1=CC1[C]=O'),
    E0 = (356.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22362,0.0644014,-7.67695e-05,5.14376e-08,-1.41422e-11,43032.8,26.907], Tmin=(100,'K'), Tmax=(879.18,'K')), NASAPolynomial(coeffs=[9.57788,0.0263906,-1.19156e-05,2.25814e-09,-1.57257e-13,41563.9,-12.3215], Tmin=(879.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(ROOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH]=C1COOC1[C]=O(20757)',
    structure = SMILES('[CH]=C1COOC1[C]=O'),
    E0 = (262.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65428,0.0350203,2.69652e-05,-6.37612e-08,2.7451e-11,31635.7,27.8246], Tmin=(100,'K'), Tmax=(975.272,'K')), NASAPolynomial(coeffs=[16.1854,0.0158691,-5.78854e-06,1.15227e-09,-8.9424e-14,26877.8,-51.7774], Tmin=(975.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[O][C]=CC1=COOC1(19206)',
    structure = SMILES('O=[C]C=C1[CH]OOC1'),
    E0 = (132.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41703,0.0418299,9.06614e-06,-4.23683e-08,1.8468e-11,16090.9,26.0153], Tmin=(100,'K'), Tmax=(1031.28,'K')), NASAPolynomial(coeffs=[16.0448,0.0196891,-9.04904e-06,1.87051e-09,-1.41618e-13,11234.1,-53.9252], Tmin=(1031.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopentane) + radical(C=CCJO) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C1COC1[C]=O(19022)',
    structure = SMILES('[CH]=C1COC1[C]=O'),
    E0 = (293.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99564,0.0303921,2.11774e-05,-5.17504e-08,2.22855e-11,35359.3,24.7271], Tmin=(100,'K'), Tmax=(978.668,'K')), NASAPolynomial(coeffs=[14.0545,0.0143921,-5.31827e-06,1.05231e-09,-8.08101e-14,31404.8,-41.3343], Tmin=(978.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'O=[C][CH]C1=COC1(20758)',
    structure = SMILES('O=[C]C=C1[CH]OC1'),
    E0 = (157.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12597,0.0412669,-2.39415e-05,5.7794e-09,-5.12132e-13,19026.2,21.9606], Tmin=(100,'K'), Tmax=(2142.18,'K')), NASAPolynomial(coeffs=[19.5542,0.0126196,-6.60998e-06,1.23461e-09,-8.0815e-14,10665.4,-77.4845], Tmin=(2142.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutane) + radical(C=CCJ(O)C) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C=C[C]=O(10456)',
    structure = SMILES('[CH]=C=C[C]=O'),
    E0 = (368.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1855,455,950,3120,650,792.5,1650,184.363,184.368,184.447,1693.07],'cm^-1')),
        HinderedRotor(inertia=(0.00837741,'amu*angstrom^2'), symmetry=1, barrier=(2.94921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0581,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58161,0.0342376,-4.2892e-05,3.11333e-08,-9.52264e-12,44355.1,15.702], Tmin=(100,'K'), Tmax=(783.956,'K')), NASAPolynomial(coeffs=[6.13378,0.0161129,-8.21217e-06,1.64147e-09,-1.17678e-13,43798.2,-0.570526], Tmin=(783.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C([CH]O[O])C[C]=O(19752)',
    structure = SMILES('[CH]=C([CH]O[O])C[C]=O'),
    E0 = (462.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0176885,'amu*angstrom^2'), symmetry=1, barrier=(3.50219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105172,'amu*angstrom^2'), symmetry=1, barrier=(20.8126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105127,'amu*angstrom^2'), symmetry=1, barrier=(20.8124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105128,'amu*angstrom^2'), symmetry=1, barrier=(20.8122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.95898,0.0629485,-5.73111e-05,2.49798e-08,-4.30733e-12,55742.3,30.436], Tmin=(100,'K'), Tmax=(1385.24,'K')), NASAPolynomial(coeffs=[16.3305,0.0185619,-9.24712e-06,1.84827e-09,-1.32661e-13,51483.7,-48.7317], Tmin=(1385.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'C=C([CH][C]=O)[CH]O[O](20759)',
    structure = SMILES('[CH2]C([CH]O[O])=C[C]=O'),
    E0 = (343.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72205,0.0767824,-9.37916e-05,5.90503e-08,-1.49814e-11,41412.1,28.4603], Tmin=(100,'K'), Tmax=(953.361,'K')), NASAPolynomial(coeffs=[12.9404,0.0255185,-1.31346e-05,2.6489e-09,-1.91433e-13,39082.3,-29.9028], Tmin=(953.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=CCJ=O)"""),
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
    label = '[CH]=C([CH2])C=C=O(16875)',
    structure = SMILES('[CH]C(=C)C=C=O'),
    E0 = (357.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.17559,'amu*angstrom^2'), symmetry=1, barrier=(50.0211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18674,'amu*angstrom^2'), symmetry=1, barrier=(50.2774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65598,0.0548559,-6.1877e-05,4.30589e-08,-1.28209e-11,43023.8,20.3781], Tmin=(100,'K'), Tmax=(802.514,'K')), NASAPolynomial(coeffs=[6.905,0.0286946,-1.29816e-05,2.443e-09,-1.69032e-13,42181.3,-3.7909], Tmin=(802.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]CO[O](20242)',
    structure = SMILES('[CH][C]CO[O]'),
    E0 = (847.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,404.609,404.615,404.632,3167.93,3167.98],'cm^-1')),
        HinderedRotor(inertia=(0.0766569,'amu*angstrom^2'), symmetry=1, barrier=(8.89488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622243,'amu*angstrom^2'), symmetry=1, barrier=(72.2622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0765029,'amu*angstrom^2'), symmetry=1, barrier=(8.89488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01219,0.0466378,-7.05672e-05,5.68391e-08,-1.80045e-11,102004,20.8821], Tmin=(100,'K'), Tmax=(805.913,'K')), NASAPolynomial(coeffs=[8.52167,0.0130199,-5.55933e-06,1.04759e-09,-7.2263e-14,100998,-8.85402], Tmin=(805.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1(C=C=O)COO1(20760)',
    structure = SMILES('[CH]C1(C=C=O)COO1'),
    E0 = (353.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.243936,0.0899555,-0.000144495,1.22285e-07,-4.0285e-11,42673.1,21.6415], Tmin=(100,'K'), Tmax=(855.15,'K')), NASAPolynomial(coeffs=[10.8261,0.0277421,-1.3065e-05,2.43672e-09,-1.64587e-13,41328.2,-25.0371], Tmin=(855.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(12dioxetane) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1(C=C=O)CO1(20761)',
    structure = SMILES('[CH]C1(C=C=O)CO1'),
    E0 = (302.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416699,0.0780005,-0.000110737,7.81948e-08,-2.05067e-11,36505.2,20.7896], Tmin=(100,'K'), Tmax=(1106.92,'K')), NASAPolynomial(coeffs=[14.1114,0.0130431,-1.74965e-06,-7.11168e-11,2.13114e-14,34421.1,-42.3899], Tmin=(1106.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(Ethylene_oxide) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C=O)=COO(20762)',
    structure = SMILES('[CH]=C([CH]OO)C=C=O'),
    E0 = (269.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586747,0.0785683,-9.94329e-05,6.46444e-08,-1.67385e-11,32487.6,28.6284], Tmin=(100,'K'), Tmax=(941.263,'K')), NASAPolynomial(coeffs=[13.6026,0.0232569,-1.12897e-05,2.21611e-09,-1.57708e-13,30037.3,-33.3776], Tmin=(941.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]([C]=C[O])CO[O](20763)',
    structure = SMILES('[CH][C]([C]=C[O])CO[O]'),
    E0 = (786.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445543,0.0750759,-8.6549e-05,4.98922e-08,-1.12115e-11,94750,32.0024], Tmin=(100,'K'), Tmax=(1093.83,'K')), NASAPolynomial(coeffs=[16.3766,0.0168169,-6.65503e-06,1.19741e-09,-8.18236e-14,91264.9,-46.284], Tmin=(1093.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(786.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(CCJ(C)CO) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1(CO[O])C=[C]O1(20764)',
    structure = SMILES('[CH]C1(CO[O])C=[C]O1'),
    E0 = (665.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.48118,0.0761664,-7.9787e-05,3.90052e-08,-7.03355e-12,80200.6,32.7616], Tmin=(100,'K'), Tmax=(1587.3,'K')), NASAPolynomial(coeffs=[22.2944,0.00489836,6.71958e-07,-2.94059e-10,2.33803e-14,74718,-82.135], Tmin=(1587.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(ROOJ) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH][C]1[CH]C(=O)OOC1(20765)',
    structure = SMILES('[CH][C]1[CH]C(=O)OOC1'),
    E0 = (376.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08215,0.0152809,9.83403e-05,-1.38527e-07,5.24877e-11,45336.1,28.0086], Tmin=(100,'K'), Tmax=(993.756,'K')), NASAPolynomial(coeffs=[18.2953,0.0179988,-8.37006e-06,1.8956e-09,-1.55304e-13,38757.1,-66.9976], Tmin=(993.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(C2CJCOOH) + radical(CCJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1C=[C]OOOC1(20766)',
    structure = SMILES('[CH]C1=C[C]OOOC1'),
    E0 = (813.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33271,0.00473469,0.000203525,-3.10982e-07,1.32587e-10,98031,27.683], Tmin=(100,'K'), Tmax=(905.335,'K')), NASAPolynomial(coeffs=[40.3265,-0.0262955,2.09008e-05,-4.16328e-09,2.71973e-13,85181.7,-188.532], Tmin=(905.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(813.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(CH2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C=O)=CO[O](20767)',
    structure = SMILES('[CH]=C([CH]O[O])C=C=O'),
    E0 = (421.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.576755,'amu*angstrom^2'), symmetry=1, barrier=(13.2607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.576938,'amu*angstrom^2'), symmetry=1, barrier=(13.2649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43479,'amu*angstrom^2'), symmetry=1, barrier=(32.9887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705411,0.0767663,-0.000110778,8.31668e-08,-2.47252e-11,50764.2,28.9562], Tmin=(100,'K'), Tmax=(825.735,'K')), NASAPolynomial(coeffs=[12.1377,0.0213866,-1.01781e-05,1.94708e-09,-1.35246e-13,48876.2,-24.0092], Tmin=(825.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C([CH]O[O])C=C=O(20768)',
    structure = SMILES('[CH]C([CH]O[O])C=C=O'),
    E0 = (565.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.368776,0.0876204,-0.000142515,1.22549e-07,-4.09754e-11,68149.6,31.0607], Tmin=(100,'K'), Tmax=(852.373,'K')), NASAPolynomial(coeffs=[10.0749,0.0281901,-1.35009e-05,2.53616e-09,-1.71994e-13,66999.3,-11.2571], Tmin=(852.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([C]=C=O)CO[O](20769)',
    structure = SMILES('[CH]C([C]=C=O)CO[O]'),
    E0 = (582.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,2120,512.5,787.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335725,0.0886103,-0.00014501,1.25178e-07,-4.19433e-11,70220.4,30.8605], Tmin=(100,'K'), Tmax=(854.216,'K')), NASAPolynomial(coeffs=[10.041,0.0284904,-1.36727e-05,2.56794e-09,-1.74022e-13,69097.6,-11.2989], Tmin=(854.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CC(C)CJ=C=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C=O)CO[O](19432)',
    structure = SMILES('[CH]=C(C=C=O)CO[O]'),
    E0 = (303.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.484939,'amu*angstrom^2'), symmetry=1, barrier=(11.1497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.483989,'amu*angstrom^2'), symmetry=1, barrier=(11.1279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.484492,'amu*angstrom^2'), symmetry=1, barrier=(11.1394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441505,0.0886202,-0.000154245,1.39363e-07,-4.79658e-11,36660.4,29.2974], Tmin=(100,'K'), Tmax=(868.944,'K')), NASAPolynomial(coeffs=[8.36739,0.0299059,-1.45179e-05,2.72182e-09,-1.83299e-13,36122.2,-2.99809], Tmin=(868.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[O]OCC1=CC([O])=C1(20770)',
    structure = SMILES('[O]OCC1=CC(=O)[CH]1'),
    E0 = (232.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52755,0.0441452,-1.09241e-05,-1.3562e-08,6.41877e-12,28059.4,23.8916], Tmin=(100,'K'), Tmax=(1190.47,'K')), NASAPolynomial(coeffs=[14.2082,0.0233278,-1.21493e-05,2.49918e-09,-1.82886e-13,23496.1,-45.981], Tmin=(1190.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(ROOJ) + radical(CCJC=O)"""),
)

species(
    label = 'C=C([C]=[C][O])CO[O](20771)',
    structure = SMILES('[CH2]C(=[C][C]=O)CO[O]'),
    E0 = (470.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.648842,'amu*angstrom^2'), symmetry=1, barrier=(14.9182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648035,'amu*angstrom^2'), symmetry=1, barrier=(14.8996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648263,'amu*angstrom^2'), symmetry=1, barrier=(14.9048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648708,'amu*angstrom^2'), symmetry=1, barrier=(14.9151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0701449,0.0966366,-0.000168848,1.52271e-07,-5.32215e-11,56664,30.8881], Tmin=(100,'K'), Tmax=(823.531,'K')), NASAPolynomial(coeffs=[10.1989,0.0301383,-1.62131e-05,3.19854e-09,-2.23063e-13,55582.4,-12.4487], Tmin=(823.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C([C]=[C]O)CO[O](20772)',
    structure = SMILES('[CH]C(=C=[C]O)CO[O]'),
    E0 = (579.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,540,610,2055,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.21835,'amu*angstrom^2'), symmetry=1, barrier=(51.0043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21851,'amu*angstrom^2'), symmetry=1, barrier=(51.0078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21844,'amu*angstrom^2'), symmetry=1, barrier=(51.0063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.21795,'amu*angstrom^2'), symmetry=1, barrier=(50.9951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271097,0.086192,-0.000113295,7.64718e-08,-1.92494e-11,69847.2,31.556], Tmin=(100,'K'), Tmax=(733.12,'K')), NASAPolynomial(coeffs=[13.0491,0.0263709,-1.11486e-05,2.00007e-09,-1.3367e-13,67707.7,-27.9377], Tmin=(733.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(AllylJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([CH]O[O])C=[C]O(20773)',
    structure = SMILES('[CH]=C([CH]O[O])C=[C]O'),
    E0 = (522.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.999061,'amu*angstrom^2'), symmetry=1, barrier=(22.9704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999124,'amu*angstrom^2'), symmetry=1, barrier=(22.9718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999776,'amu*angstrom^2'), symmetry=1, barrier=(22.9868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00013,'amu*angstrom^2'), symmetry=1, barrier=(22.995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.510912,0.0830156,-9.65793e-05,5.20234e-08,-1.03721e-11,62964.1,34.0493], Tmin=(100,'K'), Tmax=(1400.66,'K')), NASAPolynomial(coeffs=[23.4518,0.0030991,1.30356e-06,-4.19128e-10,3.29889e-14,57377.9,-85.61], Tmin=(1400.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (445.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (690.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1149.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (452.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (500.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (470.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (790.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (627.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (525.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (607.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (556.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (720.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (517.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (857.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (774.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (901.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (644.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (604.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (648.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (581.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (740.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1047.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (955.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (448.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (452.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (452.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (536.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (528.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (683.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (607.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (638.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (590.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (445.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1068.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (453.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (576.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (523.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (809.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (665.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (508.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (813.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (640.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (470.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (722.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (760.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (445.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (453.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (754.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (654.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (772.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (548.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['HCCO(2227)', 'C#CCO[O](16808)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH]C(=C[C]=O)C[O](20739)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', '[CH]C(=[CH])CO[O](20211)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C1=CC(=O)OOC1(20740)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['O(4)', '[CH]C1=CC(=O)OC1(20741)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_SDS;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C(=C=C=O)COO(20742)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C([CH][C]=O)[CH]O[O](20743)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C1(CO[O])[CH]C1=O(20744)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH][C]1COOC1[C]=O(20745)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.3498e+10,'s^-1'), n=0.261395, Ea=(79.9502,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 72.9 to 80.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C1([CH][C]=O)COO1(20746)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(162.34,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 159.6 to 162.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO(2039)', '[CH]C(=[CH])CO[O](20211)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.02e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C(=C=C=O)CO[O](20747)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cdd;HJ] for rate rule [Ca_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O2(2)', '[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.40817e+07,'m^3/(mol*s)'), n=-0.0872803, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_Ext-5R!H-R_Ext-5R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_Ext-5R!H-R_Ext-5R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_Ext-5R!H-R_Ext-5R!H-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]O[O](46)', '[CH][C]=C[C]=O(10445)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH]C([CH]O[O])=C[C]=O(20748)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]C(=[C][C]=O)CO[O](20749)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(=[C]C=O)CO[O](20750)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C([CH]OO)=C[C]=O(20751)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 336 used for R3H_SS_O;O_rad_out;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C([CH]O[O])=CC=O(20752)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.60962e+07,'s^-1'), n=1.28849, Ea=(203.441,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_1H] for rate rule [R4H_SDS;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]C(=[C][C]=O)COO(20753)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]C([C]=O)CO[O](20754)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C]CO[O](16805)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[C]=C([CH][C]=O)CO[O](20755)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[O]OCC1=CC1[C]=O(20756)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]=C1COOC1[C]=O(20757)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[O][C]=CC1=COOC1(19206)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['O(4)', '[CH]=C1COC1[C]=O(19022)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(91.0359,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation
Ea raised from 87.9 to 91.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['O(4)', 'O=[C][CH]C1=COC1(20758)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO] for rate rule [R3OO_DS;Cd_pri_rad_in;OOJ]
Euclidian distance = 2.449489742783178
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C][O](6861)', 'C#CCO[O](16808)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]O[O](46)', '[CH]=C=C[C]=O(10456)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(182.434,'m^3/(mol*s)'), n=0.88, Ea=(33.1163,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-OsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([CH]O[O])C[C]=O(19752)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00703183,'s^-1'), n=4.63833, Ea=(175.937,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeO;XH_out] + [R3H_SS_2Cd;C_rad_out_1H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['C=C([CH][C]=O)[CH]O[O](20759)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['O2(S)(666)', '[CH]=C([CH2])C=C=O(16875)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction34',
    reactants = ['HCCO(2227)', '[CH][C]CO[O](20242)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C1(C=C=O)COO1(20760)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['O(4)', '[CH]C1(C=C=O)CO1(20761)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C(C=C=O)=COO(20762)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH][C]([C]=C[O])CO[O](20763)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]C1(CO[O])C=[C]O1(20764)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.73504e+09,'s^-1'), n=0.685238, Ea=(220.064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 216.3 to 220.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH][C]1[CH]C(=O)OOC1(20765)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH][C]1C=[C]OOOC1(20766)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.34561e+10,'s^-1'), n=0.38908, Ea=(368.596,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;multiplebond_intra;radadd_intra] for rate rule [R7;multiplebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 358.2 to 368.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(3)', '[CH]C(C=C=O)=CO[O](20767)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction43',
    reactants = ['O2(2)', '[CH]=C([CH2])C=C=O(16875)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.000167406,'m^3/(mol*s)'), n=2.98833, Ea=(122.452,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O2b]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]C([CH]O[O])C=C=O(20768)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(0.0026285,'s^-1'), n=4.42, Ea=(156.9,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C([C]=C=O)CO[O](20769)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[CH]=C(C=C=O)CO[O](19432)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=C(C=[C][O])CO[O](19438)'],
    products = ['[O]OCC1=CC([O])=C1(20770)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['HCCO(2227)', '[CH]=[C]CO[O](16805)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C([C]=[C][O])CO[O](20771)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(8.97864e+06,'s^-1'), n=1.84533, Ea=(184.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SD;Y_rad_out;Cd_H_out_singleH] + [R3H_SD;Cd_rad_out;Cd_H_out_single] for rate rule [R3H_SD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=C([C]=[C]O)CO[O](20772)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C([CH]O[O])C=[C]O(20773)'],
    products = ['[CH]=C(C=[C][O])CO[O](19438)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_1H;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4459',
    isomers = [
        '[CH]=C(C=[C][O])CO[O](19438)',
    ],
    reactants = [
        ('HCCO(2227)', 'C#CCO[O](16808)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4459',
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


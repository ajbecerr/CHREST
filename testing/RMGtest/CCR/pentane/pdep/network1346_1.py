species(
    label = '[CH2]C([O])C[CH]C=C(4456)',
    structure = SMILES('[CH2]C=CCC([CH2])[O]'),
    E0 = (344.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,494.72,494.78,494.786],'cm^-1')),
        HinderedRotor(inertia=(0.0150659,'amu*angstrom^2'), symmetry=1, barrier=(2.61747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113847,'amu*angstrom^2'), symmetry=1, barrier=(2.61757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990383,'amu*angstrom^2'), symmetry=1, barrier=(22.7709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990392,'amu*angstrom^2'), symmetry=1, barrier=(22.7711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416709,0.0694634,-5.18839e-05,1.76983e-08,-1.77948e-12,41588.8,30.0755], Tmin=(100,'K'), Tmax=(1146.27,'K')), NASAPolynomial(coeffs=[15.6501,0.0271125,-1.06061e-05,1.91654e-09,-1.31405e-13,37386.5,-48.5933], Tmin=(1146.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = 'vinoxy(1351)',
    structure = SMILES('[CH2]C=O'),
    E0 = (6.53237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1127.76,1127.76,1127.77,1945],'cm^-1')),
        HinderedRotor(inertia=(0.00342653,'amu*angstrom^2'), symmetry=1, barrier=(3.09254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51396,0.0045106,2.69938e-05,-3.71579e-08,1.42465e-11,808.756,8.86322], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[6.49564,0.00770518,-2.52913e-06,4.69083e-10,-3.5128e-14,-479.661,-9.13857], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.53237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'butadiene13(1350)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30711,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80606,0.0102576,6.17291e-05,-9.01684e-08,3.59136e-11,11658.5,12.0619], Tmin=(100,'K'), Tmax=(946.037,'K')), NASAPolynomial(coeffs=[12.4692,0.0100558,-2.41229e-06,4.5713e-10,-3.93205e-14,8010.87,-43.6362], Tmin=(946.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH][CH]CCC=O(4258)',
    structure = SMILES('[CH2][CH][CH]CCC=O'),
    E0 = (319.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,180.195,1248.65,3975.57],'cm^-1')),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3572.62,'J/mol'), sigma=(6.29932,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.04 K, Pc=32.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96258,0.0421737,-1.60296e-05,1.66613e-09,4.99658e-14,38383.4,23.0896], Tmin=(100,'K'), Tmax=(2776.37,'K')), NASAPolynomial(coeffs=[48.7984,-0.00559881,6.91183e-07,-1.66324e-10,1.84404e-14,7003.1,-250.678], Tmin=(2776.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]CC=C[CH2](652)',
    structure = SMILES('[CH2][CH]CC=C[CH2]'),
    E0 = (474.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,348.938,2083.29],'cm^-1')),
        HinderedRotor(inertia=(0.107987,'amu*angstrom^2'), symmetry=1, barrier=(9.35741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0235326,'amu*angstrom^2'), symmetry=1, barrier=(56.8317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108208,'amu*angstrom^2'), symmetry=1, barrier=(9.35905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0127607,'amu*angstrom^2'), symmetry=1, barrier=(30.8466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44643,0.0505193,-2.85643e-05,7.98419e-09,-9.07072e-13,57209.5,28.1169], Tmin=(100,'K'), Tmax=(1965.55,'K')), NASAPolynomial(coeffs=[13.3496,0.0262965,-1.00793e-05,1.71471e-09,-1.09676e-13,52530.1,-37.3531], Tmin=(1965.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CC[CH][O](1848)',
    structure = SMILES('[CH2]C=CC[CH][O]'),
    E0 = (345.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,368.754,368.859,1940.62],'cm^-1')),
        HinderedRotor(inertia=(0.073507,'amu*angstrom^2'), symmetry=1, barrier=(7.09342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734922,'amu*angstrom^2'), symmetry=1, barrier=(7.09448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524237,'amu*angstrom^2'), symmetry=1, barrier=(50.5887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68117,0.0552208,-5.54951e-05,3.68237e-08,-1.10436e-11,41588.9,23.6157], Tmin=(100,'K'), Tmax=(778.818,'K')), NASAPolynomial(coeffs=[5.66422,0.0347641,-1.60961e-05,3.09862e-09,-2.18027e-13,40968.5,5.39543], Tmin=(778.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CCC([CH2])[O](5608)',
    structure = SMILES('[CH]=CCC([CH2])[O]'),
    E0 = (476.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,940.396,943.456],'cm^-1')),
        HinderedRotor(inertia=(0.557788,'amu*angstrom^2'), symmetry=1, barrier=(13.2432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.495795,'amu*angstrom^2'), symmetry=1, barrier=(13.2833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.457088,'amu*angstrom^2'), symmetry=1, barrier=(13.2683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12923,0.0583126,-5.12101e-05,2.34704e-08,-4.32559e-12,57389.8,25.8314], Tmin=(100,'K'), Tmax=(1300.12,'K')), NASAPolynomial(coeffs=[13.1951,0.0211905,-8.38098e-06,1.5088e-09,-1.0262e-13,54252.4,-35.5463], Tmin=(1300.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_P)"""),
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
    label = '[CH]C([O])CC=C[CH2](5609)',
    structure = SMILES('[CH]C([O])CC=C[CH2]'),
    E0 = (581.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.534502,0.0671261,-5.09771e-05,1.74054e-08,-1.68235e-12,70043.9,30.0226], Tmin=(100,'K'), Tmax=(1134.65,'K')), NASAPolynomial(coeffs=[15.5845,0.0250573,-9.88768e-06,1.79753e-09,-1.23832e-13,65921.3,-47.6029], Tmin=(1134.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=CCC([CH2])[O](5610)',
    structure = SMILES('[CH]C=CCC([CH2])[O]'),
    E0 = (563.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597765,0.0697836,-5.27952e-05,2.10684e-08,-3.47474e-12,67940.3,30.251], Tmin=(100,'K'), Tmax=(1408.66,'K')), NASAPolynomial(coeffs=[13.6407,0.0327472,-1.33572e-05,2.40379e-09,-1.62254e-13,64265.7,-37.1424], Tmin=(1408.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC1CO1(5611)',
    structure = SMILES('[CH2]C=CCC1CO1'),
    E0 = (92.2648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13526,0.050439,9.0072e-06,-5.43769e-08,2.78742e-11,11211.9,24.0943], Tmin=(100,'K'), Tmax=(888.95,'K')), NASAPolynomial(coeffs=[14.731,0.0241681,-5.56267e-06,7.20988e-10,-4.33925e-14,7415.56,-47.654], Tmin=(888.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.2648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1CC=CCO1(4893)',
    structure = SMILES('[CH2]C1CC=CCO1'),
    E0 = (53.8409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33852,0.0451217,1.75208e-05,-5.54487e-08,2.52593e-11,6583.64,19.9708], Tmin=(100,'K'), Tmax=(938.251,'K')), NASAPolynomial(coeffs=[13.2887,0.0278579,-8.72901e-06,1.46539e-09,-1.00833e-13,2858.6,-44.8212], Tmin=(938.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.8409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]C1CC=CCC1(4897)',
    structure = SMILES('[O]C1CC=CCC1'),
    E0 = (36.6215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80585,0.0274301,7.08305e-05,-1.08341e-07,4.22363e-11,4502.22,20.5583], Tmin=(100,'K'), Tmax=(974.041,'K')), NASAPolynomial(coeffs=[14.3294,0.0280247,-1.02009e-05,1.95391e-09,-1.46031e-13,-405.382,-52.1998], Tmin=(974.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.6215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C=CCC(=C)O(5612)',
    structure = SMILES('[CH2]C=CCC(=C)O'),
    E0 = (-10.3352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574916,0.0600404,-1.13056e-05,-3.61068e-08,2.03694e-11,-1105.84,27.2609], Tmin=(100,'K'), Tmax=(954.147,'K')), NASAPolynomial(coeffs=[19.1321,0.0194546,-5.99854e-06,1.05736e-09,-7.7201e-14,-6340.89,-70.2716], Tmin=(954.147,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.3352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC(C)=O(5613)',
    structure = SMILES('[CH2]C=CCC(C)=O'),
    E0 = (-32.5359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16031,0.0506775,-5.76552e-06,-2.18758e-08,9.51957e-12,-3801.18,25.55], Tmin=(100,'K'), Tmax=(1137.21,'K')), NASAPolynomial(coeffs=[13.505,0.0324213,-1.4878e-05,2.92468e-09,-2.10104e-13,-8236.12,-42.7479], Tmin=(1137.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.5359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([O])CC=CC(2558)',
    structure = SMILES('C=C([O])CC=CC'),
    E0 = (-24.0296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,234.646,234.758,234.823],'cm^-1')),
        HinderedRotor(inertia=(0.00306028,'amu*angstrom^2'), symmetry=1, barrier=(0.119628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251283,'amu*angstrom^2'), symmetry=1, barrier=(9.83273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251422,'amu*angstrom^2'), symmetry=1, barrier=(9.83137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885366,0.0592119,-2.96115e-05,-2.80559e-09,5.1329e-12,-2769.91,27.3679], Tmin=(100,'K'), Tmax=(1031.28,'K')), NASAPolynomial(coeffs=[13.5743,0.0280751,-1.06196e-05,1.91648e-09,-1.32745e-13,-6348.48,-38.9009], Tmin=(1031.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.0296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)CC=C=C(4896)',
    structure = SMILES('[CH2]C(O)CC=C=C'),
    E0 = (139.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580259,0.0722854,-6.4144e-05,3.03884e-08,-5.87191e-12,16890.9,28.4491], Tmin=(100,'K'), Tmax=(1231.25,'K')), NASAPolynomial(coeffs=[13.7248,0.0295826,-1.21204e-05,2.21994e-09,-1.52449e-13,13654.1,-37.7], Tmin=(1231.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO)"""),
)

species(
    label = 'C=C=CCC(C)[O](4900)',
    structure = SMILES('C=C=CCC(C)[O]'),
    E0 = (158.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600231,0.0663907,-4.90773e-05,1.86243e-08,-2.86266e-12,19152.4,27.5638], Tmin=(100,'K'), Tmax=(1527.32,'K')), NASAPolynomial(coeffs=[15.4781,0.0274258,-1.08092e-05,1.92041e-09,-1.28461e-13,14607.7,-50.5144], Tmin=(1527.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])[O](5614)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])[O]'),
    E0 = (677.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,228.938,1130.48,1456.53,2126.86],'cm^-1')),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00051,0.0692042,-6.53765e-05,3.68934e-08,-8.9716e-12,81554.8,34.6887], Tmin=(100,'K'), Tmax=(965.242,'K')), NASAPolynomial(coeffs=[8.60266,0.0377007,-1.64197e-05,3.08033e-09,-2.13975e-13,80087.2,-1.71839], Tmin=(965.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C]([CH2])[O](5615)',
    structure = SMILES('[CH2][CH]CC[C]([CH2])[O]'),
    E0 = (653.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,232.124,807.996,1204.69,1799.26],'cm^-1')),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430242,0.0859354,-0.000121676,1.03884e-07,-3.57202e-11,78771.7,33.6524], Tmin=(100,'K'), Tmax=(837.086,'K')), NASAPolynomial(coeffs=[7.05161,0.0413927,-1.87383e-05,3.48865e-09,-2.37434e-13,78115.2,5.58572], Tmin=(837.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[C]([CH2])[O](5616)',
    structure = SMILES('[CH2]C[CH]C[C]([CH2])[O]'),
    E0 = (653.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,185.833,675.106,843.071,3300.92],'cm^-1')),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467153,0.0864175,-0.000126902,1.1233e-07,-3.95982e-11,78770.5,33.3478], Tmin=(100,'K'), Tmax=(841.968,'K')), NASAPolynomial(coeffs=[5.97125,0.0433223,-1.99352e-05,3.72974e-09,-2.54111e-13,78444.3,11.3076], Tmin=(841.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CC1[CH]C1(5617)',
    structure = SMILES('[CH2]C([O])CC1[CH]C1'),
    E0 = (457.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987947,0.0550204,-1.47578e-05,-1.82047e-08,1.02782e-11,55092.5,29.3262], Tmin=(100,'K'), Tmax=(1027.59,'K')), NASAPolynomial(coeffs=[13.8437,0.028801,-1.12598e-05,2.08698e-09,-1.47363e-13,51192.7,-39.1653], Tmin=(1027.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(cyclopropane) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])O1(5566)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])O1'),
    E0 = (349.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01151,0.0539415,-2.01962e-07,-4.75318e-08,2.65035e-11,42155.7,23.3552], Tmin=(100,'K'), Tmax=(875.175,'K')), NASAPolynomial(coeffs=[15.5387,0.0219248,-4.25282e-06,4.40663e-10,-2.23713e-14,38296.2,-52.3154], Tmin=(875.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1[CH]CC([O])C1(5530)',
    structure = SMILES('[CH2]C1[CH]CC([O])C1'),
    E0 = (333.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86881,0.0289138,6.00029e-05,-9.61963e-08,3.85073e-11,40146,27.0332], Tmin=(100,'K'), Tmax=(953.699,'K')), NASAPolynomial(coeffs=[12.876,0.0277933,-9.08397e-06,1.62356e-09,-1.17455e-13,35997.9,-36.2885], Tmin=(953.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])O1(5565)',
    structure = SMILES('[CH2][CH]C1CC([CH2])O1'),
    E0 = (422.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05909,0.0536278,-2.12268e-06,-4.32976e-08,2.43253e-11,50957,27.4182], Tmin=(100,'K'), Tmax=(877.909,'K')), NASAPolynomial(coeffs=[14.7512,0.0233221,-5.15316e-06,6.25977e-10,-3.54616e-14,47316.7,-43.8959], Tmin=(877.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJCO) + radical(CJC(C)OC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([O])C1(5529)',
    structure = SMILES('[CH2][CH]C1CC([O])C1'),
    E0 = (414.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51773,0.0394184,2.88431e-05,-6.063e-08,2.43389e-11,50009.6,28.1078], Tmin=(100,'K'), Tmax=(1007.43,'K')), NASAPolynomial(coeffs=[12.879,0.0306433,-1.21923e-05,2.32672e-09,-1.68609e-13,45876.7,-35.9393], Tmin=(1007.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C=C[CH]CC(=C)[O](4317)',
    structure = SMILES('[CH2]C=CCC(=C)[O]'),
    E0 = (127.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,367.998,368.338,368.342],'cm^-1')),
        HinderedRotor(inertia=(0.00124219,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170813,'amu*angstrom^2'), symmetry=1, barrier=(16.4576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170662,'amu*angstrom^2'), symmetry=1, barrier=(16.4578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3874.43,'J/mol'), sigma=(6.44487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.18 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924688,0.0562639,-1.95361e-05,-1.80404e-08,1.18946e-11,15451.9,27.5064], Tmin=(100,'K'), Tmax=(976.244,'K')), NASAPolynomial(coeffs=[15.5,0.0226347,-7.95274e-06,1.42522e-09,-1.00663e-13,11362.9,-48.8281], Tmin=(976.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([O])CC=C=C(5618)',
    structure = SMILES('[CH2]C([O])CC=C=C'),
    E0 = (369.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,180,180,224.729],'cm^-1')),
        HinderedRotor(inertia=(0.043049,'amu*angstrom^2'), symmetry=1, barrier=(14.9271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.649231,'amu*angstrom^2'), symmetry=1, barrier=(14.9271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0430057,'amu*angstrom^2'), symmetry=1, barrier=(14.9272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700471,0.0687853,-5.99085e-05,2.7415e-08,-5.0934e-12,44593.3,27.826], Tmin=(100,'K'), Tmax=(1280,'K')), NASAPolynomial(coeffs=[13.9927,0.0272473,-1.12314e-05,2.06247e-09,-1.41759e-13,41190.4,-39.5829], Tmin=(1280,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C=CCC=O(2611)',
    structure = SMILES('[CH2]C=CCC=O'),
    E0 = (22.3348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,911.918],'cm^-1')),
        HinderedRotor(inertia=(0.0349182,'amu*angstrom^2'), symmetry=1, barrier=(20.6164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896283,'amu*angstrom^2'), symmetry=1, barrier=(20.6073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339207,'amu*angstrom^2'), symmetry=1, barrier=(20.6178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90907,0.0343719,1.38988e-05,-3.66949e-08,1.41585e-11,2771.34,21.2658], Tmin=(100,'K'), Tmax=(1086.76,'K')), NASAPolynomial(coeffs=[11.6938,0.0257614,-1.20408e-05,2.4205e-09,-1.77265e-13,-973.667,-34.1995], Tmin=(1086.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.3348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]CC=C(4889)',
    structure = SMILES('[CH2]C=CCC=C'),
    E0 = (204.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,272.006,272.35],'cm^-1')),
        HinderedRotor(inertia=(0.411509,'amu*angstrom^2'), symmetry=1, barrier=(21.5865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00227176,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409484,'amu*angstrom^2'), symmetry=1, barrier=(21.5978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56194,0.0414538,1.08814e-05,-4.06106e-08,1.75672e-11,24660.9,23.662], Tmin=(100,'K'), Tmax=(1001.11,'K')), NASAPolynomial(coeffs=[12.2191,0.0272829,-1.04542e-05,1.94474e-09,-1.38777e-13,21103.4,-34.8748], Tmin=(1001.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C=C(3743)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210055,'amu*angstrom^2'), symmetry=1, barrier=(25.2323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779472,'amu*angstrom^2'), symmetry=1, barrier=(93.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56316,0.022343,1.87062e-05,-3.93092e-08,1.63979e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.267,'K')), NASAPolynomial(coeffs=[9.83,0.0151965,-5.22268e-06,9.67646e-10,-7.07852e-14,30607.7,-26.9852], Tmin=(974.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][O](1556)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (367.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1667.93],'cm^-1')),
        HinderedRotor(inertia=(0.00517725,'amu*angstrom^2'), symmetry=1, barrier=(10.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22654,0.0212781,-3.59354e-05,3.91027e-08,-1.59281e-11,44278.5,12.2199], Tmin=(100,'K'), Tmax=(830.699,'K')), NASAPolynomial(coeffs=[2.17156,0.016018,-7.76586e-06,1.51127e-09,-1.05387e-13,44810.5,19.2612], Tmin=(830.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH]C=C(3739)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,230.296,230.449,230.657],'cm^-1')),
        HinderedRotor(inertia=(1.34138,'amu*angstrom^2'), symmetry=1, barrier=(50.5187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,44095,-3.44607], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])[O](1366)',
    structure = SMILES('[CH2]C([CH2])[O]'),
    E0 = (362.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,450.564],'cm^-1')),
        HinderedRotor(inertia=(0.000831827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.050104,'amu*angstrom^2'), symmetry=1, barrier=(7.38455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17517,0.0359905,-3.15384e-05,1.43889e-08,-2.59785e-12,43692.5,16.9615], Tmin=(100,'K'), Tmax=(1341.41,'K')), NASAPolynomial(coeffs=[10.4304,0.0113739,-4.01123e-06,7.08062e-10,-4.81153e-14,41477.8,-25.2898], Tmin=(1341.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CC(C)OJ) + radical(CJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])[O](5434)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])[O]'),
    E0 = (461.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,459.869,459.869,459.87],'cm^-1')),
        HinderedRotor(inertia=(0.598018,'amu*angstrom^2'), symmetry=1, barrier=(89.7449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179001,'amu*angstrom^2'), symmetry=1, barrier=(26.8625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178999,'amu*angstrom^2'), symmetry=1, barrier=(26.8625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179,'amu*angstrom^2'), symmetry=1, barrier=(26.8625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234313,0.0730297,-5.84454e-05,2.02624e-08,-1.74389e-12,55657.7,27.9272], Tmin=(100,'K'), Tmax=(1090.56,'K')), NASAPolynomial(coeffs=[17.3969,0.0237107,-9.35774e-06,1.71542e-09,-1.1942e-13,51103.8,-60.0762], Tmin=(1090.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]CC(=C)[O](4319)',
    structure = SMILES('[CH2][CH][CH]CC(=C)[O]'),
    E0 = (453.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1691.94,1692.09,1692.41,1692.43],'cm^-1')),
        HinderedRotor(inertia=(0.0775941,'amu*angstrom^2'), symmetry=1, barrier=(5.04212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0773371,'amu*angstrom^2'), symmetry=1, barrier=(5.03386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.078018,'amu*angstrom^2'), symmetry=1, barrier=(5.06108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779554,'amu*angstrom^2'), symmetry=1, barrier=(5.04509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2439,0.0656708,-7.98273e-05,6.56523e-08,-2.26796e-11,54615.6,33.1053], Tmin=(100,'K'), Tmax=(830.908,'K')), NASAPolynomial(coeffs=[4.97977,0.0392941,-1.70605e-05,3.13676e-09,-2.12872e-13,54284.4,17.5173], Tmin=(830.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])[O](5435)',
    structure = SMILES('[CH2]C=[C]CC([CH2])[O]'),
    E0 = (582.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,883.629,883.65,883.679],'cm^-1')),
        HinderedRotor(inertia=(0.585325,'amu*angstrom^2'), symmetry=1, barrier=(13.4578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19478,'amu*angstrom^2'), symmetry=1, barrier=(4.47837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585325,'amu*angstrom^2'), symmetry=1, barrier=(13.4578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15921,'amu*angstrom^2'), symmetry=1, barrier=(26.6525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706448,0.0703285,-6.39422e-05,3.09242e-08,-6.09694e-12,70178,29.5052], Tmin=(100,'K'), Tmax=(1207.87,'K')), NASAPolynomial(coeffs=[13.3589,0.0284273,-1.19058e-05,2.20268e-09,-1.52133e-13,67121.6,-33.9249], Tmin=(1207.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])[O](5437)',
    structure = SMILES('[CH2][C]=CCC([CH2])[O]'),
    E0 = (582.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,593.525,593.955],'cm^-1')),
        HinderedRotor(inertia=(0.104984,'amu*angstrom^2'), symmetry=1, barrier=(2.41379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104795,'amu*angstrom^2'), symmetry=1, barrier=(2.40945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0653449,'amu*angstrom^2'), symmetry=1, barrier=(16.3428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13532,'amu*angstrom^2'), symmetry=1, barrier=(26.1032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706448,0.0703285,-6.39422e-05,3.09242e-08,-6.09694e-12,70178,29.5052], Tmin=(100,'K'), Tmax=(1207.87,'K')), NASAPolynomial(coeffs=[13.3589,0.0284273,-1.19058e-05,2.20268e-09,-1.52133e-13,67121.6,-33.9249], Tmin=(1207.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH][CH]CC(=C)O(4326)',
    structure = SMILES('[CH2][CH][CH]CC(=C)O'),
    E0 = (315.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,203.845,522.076,2865.65],'cm^-1')),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03723,0.0676371,-6.46233e-05,3.7456e-08,-9.31459e-12,38051.7,32.3541], Tmin=(100,'K'), Tmax=(949.571,'K')), NASAPolynomial(coeffs=[8.4528,0.0363987,-1.52762e-05,2.81001e-09,-1.9288e-13,36643.4,-3.03783], Tmin=(949.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)=O(4331)',
    structure = SMILES('[CH2][CH][CH]CC(C)=O'),
    E0 = (293.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3139,0.0642699,-7.04167e-05,5.80813e-08,-2.14531e-11,35432.6,32.134], Tmin=(100,'K'), Tmax=(769.892,'K')), NASAPolynomial(coeffs=[3.70612,0.0452824,-2.06448e-05,3.91782e-09,-2.72106e-13,35258.6,22.4809], Tmin=(769.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([O])CC=[C]C(5619)',
    structure = SMILES('[CH2]C([O])CC=[C]C'),
    E0 = (430.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908512,0.070374,-6.36872e-05,3.26535e-08,-7.10365e-12,51945.9,28.5055], Tmin=(100,'K'), Tmax=(1075.57,'K')), NASAPolynomial(coeffs=[10.1358,0.036058,-1.58293e-05,2.98973e-09,-2.08685e-13,49961,-16.6828], Tmin=(1075.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])O(5620)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])O'),
    E0 = (231.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0720449,0.0786114,-6.94006e-05,3.11102e-08,-5.51262e-12,27963.5,29.2252], Tmin=(100,'K'), Tmax=(1366.72,'K')), NASAPolynomial(coeffs=[18.8038,0.0233672,-8.76926e-06,1.53519e-09,-1.02792e-13,22803.9,-67.7372], Tmin=(1366.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)[O](5621)',
    structure = SMILES('[CH2]C=C[CH]C(C)[O]'),
    E0 = (249.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420245,0.0674528,-3.73873e-05,-5.16328e-10,5.06711e-12,30204.1,26.6255], Tmin=(100,'K'), Tmax=(1048.99,'K')), NASAPolynomial(coeffs=[16.3814,0.0278217,-1.10777e-05,2.05836e-09,-1.45074e-13,25687.3,-56.7091], Tmin=(1048.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([O])C[C]=CC(5622)',
    structure = SMILES('[CH2]C([O])C[C]=CC'),
    E0 = (430.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,211.083,211.158,3552],'cm^-1')),
        HinderedRotor(inertia=(0.00376494,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.533059,'amu*angstrom^2'), symmetry=1, barrier=(16.9504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538182,'amu*angstrom^2'), symmetry=1, barrier=(16.9623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534462,'amu*angstrom^2'), symmetry=1, barrier=(16.9386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908512,0.070374,-6.36872e-05,3.26535e-08,-7.10365e-12,51945.9,28.5055], Tmin=(100,'K'), Tmax=(1075.57,'K')), NASAPolynomial(coeffs=[10.1358,0.036058,-1.58293e-05,2.98973e-09,-2.08685e-13,49961,-16.6828], Tmin=(1075.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])O(5623)',
    structure = SMILES('[CH2]C=[C]CC([CH2])O'),
    E0 = (352.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.765838,'amu*angstrom^2'), symmetry=1, barrier=(17.6081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00656176,'amu*angstrom^2'), symmetry=1, barrier=(2.96051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765967,'amu*angstrom^2'), symmetry=1, barrier=(17.6111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766383,'amu*angstrom^2'), symmetry=1, barrier=(17.6207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83218,'amu*angstrom^2'), symmetry=1, barrier=(42.1255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.572637,0.0739736,-6.86168e-05,3.43785e-08,-7.04663e-12,42476.3,30.1781], Tmin=(100,'K'), Tmax=(1163.17,'K')), NASAPolynomial(coeffs=[13.1753,0.0306342,-1.27268e-05,2.34513e-09,-1.6164e-13,39544.5,-32.5272], Tmin=(1163.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(C)[O](5624)',
    structure = SMILES('[CH2]C=[C]CC(C)[O]'),
    E0 = (370.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,285.022,285.182,1341.82],'cm^-1')),
        HinderedRotor(inertia=(0.00207577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177811,'amu*angstrom^2'), symmetry=1, barrier=(10.258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357656,'amu*angstrom^2'), symmetry=1, barrier=(20.6544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357444,'amu*angstrom^2'), symmetry=1, barrier=(20.6536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666022,0.0673342,-5.14422e-05,2.04652e-08,-3.32533e-12,44734.2,29.0204], Tmin=(100,'K'), Tmax=(1440.85,'K')), NASAPolynomial(coeffs=[14.3753,0.0292745,-1.18191e-05,2.13158e-09,-1.44221e-13,40783.7,-42.1257], Tmin=(1440.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C=C[CH]C(4455)',
    structure = SMILES('[CH2]C([O])[CH]C=CC'),
    E0 = (310.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,265.529,266.44,266.545],'cm^-1')),
        HinderedRotor(inertia=(0.00237615,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488918,'amu*angstrom^2'), symmetry=1, barrier=(24.6072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.493759,'amu*angstrom^2'), symmetry=1, barrier=(24.6065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.36947,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354956,0.0740899,-6.18609e-05,2.66954e-08,-4.6568e-12,37428.8,27.2142], Tmin=(100,'K'), Tmax=(1360.75,'K')), NASAPolynomial(coeffs=[15.8174,0.0286372,-1.17569e-05,2.14813e-09,-1.46926e-13,33220.7,-52.1461], Tmin=(1360.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])O(5625)',
    structure = SMILES('[CH2][C]=CCC([CH2])O'),
    E0 = (352.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0390436,'amu*angstrom^2'), symmetry=1, barrier=(17.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000260745,'amu*angstrom^2'), symmetry=1, barrier=(2.96052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83219,'amu*angstrom^2'), symmetry=1, barrier=(42.1258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766205,'amu*angstrom^2'), symmetry=1, barrier=(17.6166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765823,'amu*angstrom^2'), symmetry=1, barrier=(17.6078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.572637,0.0739736,-6.86168e-05,3.43785e-08,-7.04663e-12,42476.3,30.1781], Tmin=(100,'K'), Tmax=(1163.17,'K')), NASAPolynomial(coeffs=[13.1753,0.0306342,-1.27268e-05,2.34513e-09,-1.6164e-13,39544.5,-32.5272], Tmin=(1163.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC(C)[O](5626)',
    structure = SMILES('[CH2][C]=CCC(C)[O]'),
    E0 = (370.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,285.022,285.182,1341.82],'cm^-1')),
        HinderedRotor(inertia=(0.00207577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177811,'amu*angstrom^2'), symmetry=1, barrier=(10.258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357656,'amu*angstrom^2'), symmetry=1, barrier=(20.6544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357444,'amu*angstrom^2'), symmetry=1, barrier=(20.6536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.666022,0.0673342,-5.14422e-05,2.04652e-08,-3.32533e-12,44734.2,29.0204], Tmin=(100,'K'), Tmax=(1440.85,'K')), NASAPolynomial(coeffs=[14.3753,0.0292745,-1.18191e-05,2.13158e-09,-1.44221e-13,40783.7,-42.1257], Tmin=(1440.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C[CH][CH]C(2560)',
    structure = SMILES('C=C([O])C[CH][CH]C'),
    E0 = (248.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,222.408,1946.36,1946.36,1946.36],'cm^-1')),
        HinderedRotor(inertia=(0.20821,'amu*angstrom^2'), symmetry=1, barrier=(7.30844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00340802,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208215,'amu*angstrom^2'), symmetry=1, barrier=(7.30845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00340791,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02514,0.0519567,-1.23522e-06,-9.02755e-08,8.70563e-11,29898,28.692], Tmin=(100,'K'), Tmax=(467.454,'K')), NASAPolynomial(coeffs=[4.96563,0.0420012,-1.80839e-05,3.34267e-09,-2.2893e-13,29457,14.9649], Tmin=(467.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([O])C([CH2])C=C(4445)',
    structure = SMILES('[CH2]C([O])C([CH2])C=C'),
    E0 = (402.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,604.318,3113.3],'cm^-1')),
        HinderedRotor(inertia=(0.0720765,'amu*angstrom^2'), symmetry=1, barrier=(18.6837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.52797,'amu*angstrom^2'), symmetry=1, barrier=(81.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721023,'amu*angstrom^2'), symmetry=1, barrier=(18.6838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139872,'amu*angstrom^2'), symmetry=1, barrier=(3.21594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3818.2,'J/mol'), sigma=(6.62498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.39 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491081,0.0696401,-5.37465e-05,1.82912e-08,-1.21122e-12,48543.4,31.4972], Tmin=(100,'K'), Tmax=(1021.56,'K')), NASAPolynomial(coeffs=[14.6069,0.0270162,-9.73121e-06,1.68648e-09,-1.13564e-13,44999.4,-40.1346], Tmin=(1021.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = 'C2H3(60)',
    structure = SMILES('[CH]=C'),
    E0 = (286.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,728.586,728.586,3521.64],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8307,-0.00147628,3.09002e-05,-3.80476e-08,1.43171e-11,34524.2,5.61949], Tmin=(100,'K'), Tmax=(933.662,'K')), NASAPolynomial(coeffs=[5.36086,0.00527345,-1.3196e-06,2.21564e-10,-1.68768e-14,33658.6,-4.76326], Tmin=(933.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]CC([CH2])[O](1572)',
    structure = SMILES('[CH]CC([CH2])[O]'),
    E0 = (573.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,883.526,883.982,884.968,4000],'cm^-1')),
        HinderedRotor(inertia=(0.209558,'amu*angstrom^2'), symmetry=1, barrier=(4.81816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732953,'amu*angstrom^2'), symmetry=1, barrier=(40.7383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625348,'amu*angstrom^2'), symmetry=1, barrier=(14.378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58742,0.0496998,-4.58357e-05,2.21998e-08,-4.31781e-12,69052.6,22.0679], Tmin=(100,'K'), Tmax=(1237.45,'K')), NASAPolynomial(coeffs=[11.5128,0.0176163,-6.94489e-06,1.24755e-09,-8.48368e-14,66596.2,-27.9307], Tmin=(1237.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC(C=C)O1(4891)',
    structure = SMILES('[CH2]C1CC(C=C)O1'),
    E0 = (145.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16246,0.0473908,2.05242e-05,-6.88981e-08,3.36092e-11,17557,24.1456], Tmin=(100,'K'), Tmax=(894.747,'K')), NASAPolynomial(coeffs=[16.1574,0.0212079,-4.06873e-06,4.55017e-10,-2.67275e-14,13238.4,-55.6668], Tmin=(894.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CC1CC([O])C1(5627)',
    structure = SMILES('C=CC1CC([O])C1'),
    E0 = (144.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60672,0.0317547,6.19162e-05,-1.02811e-07,4.13262e-11,17425.8,23.9527], Tmin=(100,'K'), Tmax=(966.216,'K')), NASAPolynomial(coeffs=[15.7607,0.0251976,-8.69122e-06,1.6478e-09,-1.24065e-13,12261.5,-56.4157], Tmin=(966.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C=CC=C(5628)',
    structure = SMILES('[CH2]C(O)C=CC=C'),
    E0 = (74.0325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491856,0.06469,-3.02503e-05,-1.30823e-08,1.13366e-11,9041.54,29.2054], Tmin=(100,'K'), Tmax=(973.469,'K')), NASAPolynomial(coeffs=[18.0091,0.0212502,-7.28966e-06,1.30908e-09,-9.34249e-14,4278.81,-61.7796], Tmin=(973.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.0325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CJCO)"""),
)

species(
    label = 'C=CC=CC(C)[O](5629)',
    structure = SMILES('C=CC=CC(C)[O]'),
    E0 = (92.8042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80615,0.0555101,-4.57115e-06,-3.73638e-08,1.91635e-11,11290,27.252], Tmin=(100,'K'), Tmax=(976.412,'K')), NASAPolynomial(coeffs=[17.1709,0.0231712,-8.20084e-06,1.51293e-09,-1.09878e-13,6440.02,-59.7785], Tmin=(976.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.8042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(=O)CCC=C(5630)',
    structure = SMILES('C=CCCC(=C)[O]'),
    E0 = (-13.0538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847324,0.0590066,-2.52037e-05,-9.83358e-09,8.18682e-12,-1447.42,27.9152], Tmin=(100,'K'), Tmax=(1004.22,'K')), NASAPolynomial(coeffs=[14.4321,0.0268413,-9.93816e-06,1.79346e-09,-1.25183e-13,-5282.38,-43.1899], Tmin=(1004.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.0538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1C[CH][CH]CO1(5579)',
    structure = SMILES('[CH2]C1C[CH][CH]CO1'),
    E0 = (324.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6963,0.0412548,1.8627e-05,-5.20293e-08,2.40383e-11,39094.2,23.797], Tmin=(100,'K'), Tmax=(890.159,'K')), NASAPolynomial(coeffs=[9.42275,0.0322363,-9.48453e-06,1.45923e-09,-9.31966e-14,36700.4,-18.2991], Tmin=(890.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(RCCJCC) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]C1C[CH][CH]CC1(5549)',
    structure = SMILES('[O]C1C[CH][CH]CC1'),
    E0 = (280.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14716,0.0326104,2.78456e-05,-4.06251e-08,1.26386e-11,33779.8,22.8964], Tmin=(100,'K'), Tmax=(1206.7,'K')), NASAPolynomial(coeffs=[5.9549,0.0458393,-2.07329e-05,3.96633e-09,-2.77812e-13,30978.8,-3.98771], Tmin=(1206.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclohexane) + radical(CC(C)OJ) + radical(cyclohexane) + radical(cyclohexane)"""),
)

species(
    label = '[CH2]C([O])C=CC=C(5631)',
    structure = SMILES('[CH2]C([O])C=CC=C'),
    E0 = (304.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,292.082,292.488,292.586],'cm^-1')),
        HinderedRotor(inertia=(0.284367,'amu*angstrom^2'), symmetry=1, barrier=(17.1699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283784,'amu*angstrom^2'), symmetry=1, barrier=(17.1743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283297,'amu*angstrom^2'), symmetry=1, barrier=(17.1731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.653667,0.060707,-2.43779e-05,-1.8079e-08,1.29277e-11,36742.1,28.4328], Tmin=(100,'K'), Tmax=(973.812,'K')), NASAPolynomial(coeffs=[17.9463,0.0194542,-6.70237e-06,1.22131e-09,-8.84206e-14,31962.2,-61.7847], Tmin=(973.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])[CH]CC=C(5587)',
    structure = SMILES('[CH2]C([O])[CH]CC=C'),
    E0 = (405.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,707.769,707.785,707.79],'cm^-1')),
        HinderedRotor(inertia=(0.00875412,'amu*angstrom^2'), symmetry=1, barrier=(3.11237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574976,'amu*angstrom^2'), symmetry=1, barrier=(13.2198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135344,'amu*angstrom^2'), symmetry=1, barrier=(3.11182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0371909,'amu*angstrom^2'), symmetry=1, barrier=(13.2198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56025,0.0665876,-4.74939e-05,1.46262e-08,-9.51909e-13,48877.7,32.214], Tmin=(100,'K'), Tmax=(1132.67,'K')), NASAPolynomial(coeffs=[14.965,0.0272444,-1.06567e-05,1.92923e-09,-1.32522e-13,44875.1,-42.3391], Tmin=(1132.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CC[C]=C(5632)',
    structure = SMILES('[CH2]C([O])CC[C]=C'),
    E0 = (443.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,316.399,316.593,1202.23,3964.94],'cm^-1')),
        HinderedRotor(inertia=(0.140197,'amu*angstrom^2'), symmetry=1, barrier=(10.0129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306529,'amu*angstrom^2'), symmetry=1, barrier=(21.9033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306067,'amu*angstrom^2'), symmetry=1, barrier=(21.8651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203,'amu*angstrom^2'), symmetry=1, barrier=(14.4032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772041,0.0712644,-6.31579e-05,3.06773e-08,-6.20726e-12,53425.3,29.4291], Tmin=(100,'K'), Tmax=(1162.35,'K')), NASAPolynomial(coeffs=[11.7624,0.0334431,-1.435e-05,2.6835e-09,-1.86311e-13,50870.4,-25.2467], Tmin=(1162.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]CCC(=C)[O](4322)',
    structure = SMILES('[CH2][CH]CCC(=C)[O]'),
    E0 = (258.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,283.133,283.138,1104.16,2570.2],'cm^-1')),
        HinderedRotor(inertia=(0.17237,'amu*angstrom^2'), symmetry=1, barrier=(9.80502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0021029,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210292,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.075,0.0641934,-5.22537e-05,2.4342e-08,-4.83332e-12,31238.6,31.147], Tmin=(100,'K'), Tmax=(1169.1,'K')), NASAPolynomial(coeffs=[9.778,0.0344163,-1.40481e-05,2.55539e-09,-1.74434e-13,29203.7,-12.1996], Tmin=(1169.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CCCC([CH2])[O](5633)',
    structure = SMILES('[CH]=CCCC([CH2])[O]'),
    E0 = (452.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542311,0.0724259,-6.28575e-05,2.88698e-08,-5.40117e-12,54549.8,30.1814], Tmin=(100,'K'), Tmax=(1269.27,'K')), NASAPolynomial(coeffs=[14.165,0.0294953,-1.2123e-05,2.22223e-09,-1.52575e-13,51091.6,-38.7883], Tmin=(1269.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]CC([CH2])O(5634)',
    structure = SMILES('[CH]C=CCC([CH2])O'),
    E0 = (333.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.495986,0.0730913,-5.64673e-05,2.34501e-08,-4.05282e-12,40237.1,30.8062], Tmin=(100,'K'), Tmax=(1344.67,'K')), NASAPolynomial(coeffs=[13.172,0.035384,-1.44045e-05,2.59605e-09,-1.75671e-13,36828.1,-34.1023], Tmin=(1344.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]CC(C)[O](5635)',
    structure = SMILES('[CH]C=CCC(C)[O]'),
    E0 = (352.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.403171,0.0683048,-4.43974e-05,1.45807e-08,-1.9468e-12,42504.1,30.3425], Tmin=(100,'K'), Tmax=(1726.55,'K')), NASAPolynomial(coeffs=[16.5164,0.0309744,-1.19654e-05,2.05787e-09,-1.33541e-13,36940,-56.1946], Tmin=(1726.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(AllylJ2_triplet)"""),
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
    E0 = (344.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (501.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (994.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (782.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (914.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (793.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (775.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (350.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (352.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (352.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (367.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (367.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (369.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (369.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (369.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (699.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (717.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (678.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (570.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (402.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (418.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (466.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (414.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (360.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (593.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (460.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (447.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (344.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (642.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (739.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (673.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (665.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (794.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (794.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (458.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (486.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (552.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (486.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (463.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (592.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (453.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (415.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (612.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (396.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (410.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (432.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (562.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (915.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (352.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (352.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (422.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (422.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (433.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (405.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (432.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (517.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (469.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (530.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (643.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (537.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (612.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (473.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (484.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['vinoxy(1351)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2][CH][CH]CCC=O(4258)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH2]C=CC[CH][O](1848)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH]=CCC([CH2])[O](5608)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C([O])CC=C[CH2](5609)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH]C=CCC([CH2])[O](5610)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C=CCC1CO1(5611)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C1CC=CCO1(4893)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[O]C1CC=CCC1(4897)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_2H] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C=CCC(=C)O(5612)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C=CCC(C)=O(5613)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['C=C([O])CC=CC(2558)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C(O)CC=C=C(4896)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['C=C=CCC(C)[O](4900)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C[CH]C([CH2])[O](5614)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]CC[C]([CH2])[O](5615)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C[CH]C[C]([CH2])[O](5616)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C([O])CC1[CH]C1(5617)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C1[CH]CC([CH2])O1(5566)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C1[CH]CC([O])C1(5530)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.47079e+07,'s^-1'), n=0.909323, Ea=(74.2834,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2][CH]C1CC([CH2])O1(5565)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2][CH]C1CC([O])C1(5529)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(70.2989,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 66.4 to 70.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C=C[CH]CC(=C)[O](4317)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]C([O])CC=C=C(5618)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(T)(20)', '[CH2]C=CCC=O(2611)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)', 'C=C[CH]CC=C(4889)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['vinoxy(1351)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.000197168,'m^3/(mol*s)'), n=2.87833, Ea=(63.4048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;CsJ-CdHH] + [CO_O;CsJ-OneDeHH] for rate rule [CO_O;CsJ-CdHH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 59.8 to 63.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.7349e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C=C(3739)', '[CH2]C([CH2])[O](1366)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.7349e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]C=C[CH]C([CH2])[O](5434)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH][CH]CC(=C)[O](4319)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2]C=[C]CC([CH2])[O](5435)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2][C]=CCC([CH2])[O](5437)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2][CH][CH]CC(=C)O(4326)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2][CH][CH]CC(C)=O(4331)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C([O])CC=[C]C(5619)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C=C[CH]C([CH2])O(5620)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.07519e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C=C[CH]C(C)[O](5621)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 85 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([O])C[C]=CC(5622)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C=[C]CC([CH2])O(5623)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C=[C]CC(C)[O](5624)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C([O])C=C[CH]C(4455)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(7.15542e+07,'s^-1'), n=1, Ea=(267.985,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] + [R4H_SDS;C_rad_out_2H;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][C]=CCC([CH2])O(5625)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C]=CCC(C)[O](5626)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['C=C([O])C[CH][CH]C(2560)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SMSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([O])C([CH2])C=C(4445)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C2H3(60)', '[CH]CC([CH2])[O](1572)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C1CC(C=C)O1(4891)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['C=CC1CC([O])C1(5627)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] + [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C(O)C=CC=C(5628)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['C=CC=CC(C)[O](5629)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C(=O)CCC=C(5630)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2]C1C[CH][CH]CO1(5579)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[O]C1C[CH][CH]CC1(5549)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['H(3)', '[CH2]C([O])C=CC=C(5631)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2][CH][O](1556)', 'butadiene13(1350)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C([O])[CH]CC=C(5587)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C([O])CC[C]=C(5632)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH2][CH]CCC(=C)[O](4322)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.32587e+10,'s^-1'), n=0.723333, Ea=(193.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH]=CCCC([CH2])[O](5633)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C([O])C[CH]C=C(4456)'],
    products = ['[CH]=C[CH]CC([CH2])O(5634)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]=C[CH]CC(C)[O](5635)'],
    products = ['[CH2]C([O])C[CH]C=C(4456)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1346',
    isomers = [
        '[CH2]C([O])C[CH]C=C(4456)',
    ],
    reactants = [
        ('vinoxy(1351)', 'butadiene13(1350)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1346',
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


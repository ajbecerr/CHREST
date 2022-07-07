species(
    label = '[CH2]C([CH][C]=C)C=C(15833)',
    structure = SMILES('[CH2][C]=CC([CH2])C=C'),
    E0 = (613.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,269.17,3233.31],'cm^-1')),
        HinderedRotor(inertia=(0.259691,'amu*angstrom^2'), symmetry=1, barrier=(13.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5476,'amu*angstrom^2'), symmetry=1, barrier=(79.5686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259692,'amu*angstrom^2'), symmetry=1, barrier=(13.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54761,'amu*angstrom^2'), symmetry=1, barrier=(79.5686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882232,0.0689034,-6.09384e-05,3.11838e-08,-6.748e-12,73908.3,27.9463], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[10.0801,0.0350785,-1.42922e-05,2.59388e-09,-1.76871e-13,71907.4,-17.2016], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56316,0.022343,1.87062e-05,-3.93092e-08,1.63979e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.267,'K')), NASAPolynomial(coeffs=[9.83,0.0151965,-5.22268e-06,9.67646e-10,-7.07852e-14,30607.7,-26.9852], Tmin=(974.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=C[CH]CC=C(17265)',
    structure = SMILES('[CH2][C]=C[CH]CC=C'),
    E0 = (558.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96128,0.0570789,-2.37796e-05,-5.6192e-09,4.99749e-12,67244.9,27.4067], Tmin=(100,'K'), Tmax=(1097.91,'K')), NASAPolynomial(coeffs=[12.9066,0.03127,-1.27165e-05,2.35633e-09,-1.64288e-13,63554.5,-36.1998], Tmin=(1097.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[CH][C]=C(15872)',
    structure = SMILES('[CH2][C]=CCC=C[CH2]'),
    E0 = (557.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,873.338],'cm^-1')),
        HinderedRotor(inertia=(0.246571,'amu*angstrom^2'), symmetry=1, barrier=(5.66916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0381156,'amu*angstrom^2'), symmetry=1, barrier=(19.8377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39198,'amu*angstrom^2'), symmetry=1, barrier=(32.0044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863749,'amu*angstrom^2'), symmetry=1, barrier=(19.8593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827604,0.0597375,-2.78019e-05,-4.2108e-09,5.17522e-12,67179.2,28.8603], Tmin=(100,'K'), Tmax=(1067.75,'K')), NASAPolynomial(coeffs=[13.8178,0.0297253,-1.18422e-05,2.18422e-09,-1.52521e-13,63341.9,-39.6401], Tmin=(1067.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH2][C]=CC=C[CH2](17118)',
    structure = SMILES('[CH2][C]=CC=C[CH2]'),
    E0 = (495.823,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,486.845],'cm^-1')),
        HinderedRotor(inertia=(0.358302,'amu*angstrom^2'), symmetry=1, barrier=(62.3806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369559,'amu*angstrom^2'), symmetry=1, barrier=(62.3686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367487,'amu*angstrom^2'), symmetry=1, barrier=(62.3799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88436,0.0352549,1.6638e-05,-4.67215e-08,2.0775e-11,59720.2,22.7901], Tmin=(100,'K'), Tmax=(947.888,'K')), NASAPolynomial(coeffs=[11.541,0.0224612,-7.35629e-06,1.26879e-09,-8.84677e-14,56633.6,-29.9054], Tmin=(947.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.823,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[CH]C(C=C)C=[C][CH2](17266)',
    structure = SMILES('[CH]C(C=C)C=[C][CH2]'),
    E0 = (856.711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.923767,0.0689395,-6.32269e-05,3.21627e-08,-6.84524e-12,103148,27.0208], Tmin=(100,'K'), Tmax=(1106.74,'K')), NASAPolynomial(coeffs=[10.9132,0.0328358,-1.42945e-05,2.68731e-09,-1.87102e-13,100937,-22.1853], Tmin=(1106.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(856.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]=CC([CH2])C=C(17267)',
    structure = SMILES('[CH][C]=CC([CH2])C=C'),
    E0 = (832.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.923291,0.0707368,-6.66183e-05,4.02263e-08,-1.07258e-11,100266,28.6345], Tmin=(100,'K'), Tmax=(879.702,'K')), NASAPolynomial(coeffs=[7.21221,0.0421403,-1.78563e-05,3.27186e-09,-2.23474e-13,99159.9,-0.89965], Tmin=(879.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(832.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC(C=C)C1(17268)',
    structure = SMILES('C=CC1[CH]C(=C)C1'),
    E0 = (329.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80872,0.0301151,5.69134e-05,-8.9833e-08,3.46825e-11,39709.8,21.7427], Tmin=(100,'K'), Tmax=(988.976,'K')), NASAPolynomial(coeffs=[12.7373,0.0306725,-1.18191e-05,2.26189e-09,-1.65992e-13,35359.3,-41.9268], Tmin=(988.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=CC(=C)C=C(17269)',
    structure = SMILES('[CH2]C(C=C)=CC=C'),
    E0 = (226.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.833978,0.0523789,1.27347e-05,-6.14115e-08,2.96777e-11,27330.4,24.1529], Tmin=(100,'K'), Tmax=(938.522,'K')), NASAPolynomial(coeffs=[18.5312,0.0203993,-5.59175e-06,9.30893e-10,-6.76453e-14,22095.1,-70.2969], Tmin=(938.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=C[C]([CH2])C[CH2](17270)',
    structure = SMILES('[CH2][C]=C[C]([CH2])C[CH2]'),
    E0 = (823.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,180,1991.87],'cm^-1')),
        HinderedRotor(inertia=(0.00553031,'amu*angstrom^2'), symmetry=1, barrier=(15.5866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0010207,'amu*angstrom^2'), symmetry=1, barrier=(2.87644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.3151,'amu*angstrom^2'), symmetry=1, barrier=(76.2206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940194,'amu*angstrom^2'), symmetry=1, barrier=(76.2292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201175,'amu*angstrom^2'), symmetry=1, barrier=(76.2048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799493,0.0643575,-4.8499e-05,2.02736e-08,-3.52284e-12,99172.9,30.8862], Tmin=(100,'K'), Tmax=(1349.31,'K')), NASAPolynomial(coeffs=[12.0531,0.0309965,-1.14125e-05,1.94997e-09,-1.2786e-13,96135.9,-26.7775], Tmin=(1349.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(823.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_T) + radical(RCCJ) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])C[CH2](17271)',
    structure = SMILES('[CH2][C]=[C]C([CH2])C[CH2]'),
    E0 = (929.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1670,1700,300,440,180,694.504],'cm^-1')),
        HinderedRotor(inertia=(0.644112,'amu*angstrom^2'), symmetry=1, barrier=(14.8094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0433058,'amu*angstrom^2'), symmetry=1, barrier=(14.8171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00109206,'amu*angstrom^2'), symmetry=1, barrier=(4.50087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0183951,'amu*angstrom^2'), symmetry=1, barrier=(75.7661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.29502,'amu*angstrom^2'), symmetry=1, barrier=(75.759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817085,0.070473,-6.52146e-05,3.49426e-08,-7.8632e-12,111899,32.3527], Tmin=(100,'K'), Tmax=(1053.53,'K')), NASAPolynomial(coeffs=[10.2918,0.0344997,-1.39965e-05,2.5322e-09,-1.72317e-13,109903,-13.8516], Tmin=(1053.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(929.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])[CH]C(17272)',
    structure = SMILES('[CH2][C]=[C]C([CH2])[CH]C'),
    E0 = (918.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,1423.6,1423.6],'cm^-1')),
        HinderedRotor(inertia=(0.00213093,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00213089,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149756,'amu*angstrom^2'), symmetry=1, barrier=(8.4072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149755,'amu*angstrom^2'), symmetry=1, barrier=(8.40721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16273,'amu*angstrom^2'), symmetry=1, barrier=(65.2743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07283,0.0672321,-6.12879e-05,3.4125e-08,-8.28106e-12,110600,32.2106], Tmin=(100,'K'), Tmax=(962.764,'K')), NASAPolynomial(coeffs=[8.0075,0.0384207,-1.63997e-05,3.04233e-09,-2.09899e-13,109265,-0.982111], Tmin=(962.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(918.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C1[CH]CC1(16131)',
    structure = SMILES('[CH2][C]=CC1[CH]CC1'),
    E0 = (627.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6491,0.036512,3.53404e-05,-6.59814e-08,2.5983e-11,75542,27.2365], Tmin=(100,'K'), Tmax=(1002.73,'K')), NASAPolynomial(coeffs=[12.0637,0.0317214,-1.24749e-05,2.3636e-09,-1.70632e-13,71605.6,-32.2497], Tmin=(1002.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC([CH2])[CH]C1(17273)',
    structure = SMILES('[CH2]C1=CC([CH2])[CH]C1'),
    E0 = (494.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63674,0.0352798,4.38829e-05,-8.00512e-08,3.27381e-11,59607.6,25.6197], Tmin=(100,'K'), Tmax=(962.178,'K')), NASAPolynomial(coeffs=[13.2686,0.0282871,-9.70047e-06,1.75489e-09,-1.26357e-13,55454.5,-39.9987], Tmin=(962.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-4) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]CC[C]=C1(17274)',
    structure = SMILES('[CH2]C1[CH]CC[C]=C1'),
    E0 = (584.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88109,0.0312833,4.70963e-05,-7.68696e-08,2.98608e-11,70382.4,24.9329], Tmin=(100,'K'), Tmax=(983.552,'K')), NASAPolynomial(coeffs=[10.9006,0.0323358,-1.20562e-05,2.23131e-09,-1.59805e-13,66783.1,-27.7097], Tmin=(983.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cs_S) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC1CC1[CH2](17275)',
    structure = SMILES('[CH2][C]=CC1CC1[CH2]'),
    E0 = (639.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27178,0.0436032,2.65995e-05,-6.85515e-08,3.05102e-11,77075.2,27.0365], Tmin=(100,'K'), Tmax=(945.252,'K')), NASAPolynomial(coeffs=[15.5964,0.0239853,-7.33021e-06,1.26444e-09,-9.06376e-14,72535.5,-50.9536], Tmin=(945.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC([CH2])C1[CH2](17276)',
    structure = SMILES('[CH2]C1[CH]C(=C)C1[CH2]'),
    E0 = (602.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59061,0.0382742,3.40594e-05,-6.98321e-08,2.9508e-11,72590.3,24.6407], Tmin=(100,'K'), Tmax=(947.882,'K')), NASAPolynomial(coeffs=[12.4593,0.0291038,-9.49723e-06,1.64319e-09,-1.14987e-13,68881.4,-35.9084], Tmin=(947.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C=[C]CC1[CH2](17277)',
    structure = SMILES('[CH2]C1C=[C]CC1[CH2]'),
    E0 = (625.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62317,0.0362555,4.1227e-05,-8.02307e-08,3.4213e-11,75281.8,25.811], Tmin=(100,'K'), Tmax=(933.785,'K')), NASAPolynomial(coeffs=[13.5096,0.0259964,-7.60484e-06,1.26088e-09,-8.82892e-14,71289.3,-40.2113], Tmin=(933.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(625.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
)

species(
    label = '[CH2][C]=CC(=C)C=C(17278)',
    structure = SMILES('[CH2]C(C=C)=C[C]=C'),
    E0 = (425.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67391,'amu*angstrom^2'), symmetry=1, barrier=(38.4864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67811,'amu*angstrom^2'), symmetry=1, barrier=(38.5831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67848,'amu*angstrom^2'), symmetry=1, barrier=(38.5915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762838,0.0588316,-1.59899e-05,-2.90974e-08,1.81929e-11,51262.2,24.0056], Tmin=(100,'K'), Tmax=(923.691,'K')), NASAPolynomial(coeffs=[17.0045,0.0204652,-5.5988e-06,8.70605e-10,-5.88226e-14,46898,-60.4441], Tmin=(923.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C#CC([CH2])C=C(17279)',
    structure = SMILES('[CH2]C#CC([CH2])C=C'),
    E0 = (516.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,2100,2250,500,550,624.83,626.681],'cm^-1')),
        HinderedRotor(inertia=(0.094472,'amu*angstrom^2'), symmetry=1, barrier=(2.17408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.612626,'amu*angstrom^2'), symmetry=1, barrier=(14.0883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00005,'amu*angstrom^2'), symmetry=1, barrier=(90.1614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.92657,'amu*angstrom^2'), symmetry=1, barrier=(90.2795,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870553,0.0650072,-5.54218e-05,2.61074e-08,-5.0587e-12,62181,27.2343], Tmin=(100,'K'), Tmax=(1226.93,'K')), NASAPolynomial(coeffs=[12.0455,0.0285749,-1.08808e-05,1.90547e-09,-1.27274e-13,59438.8,-28.9637], Tmin=(1226.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Isobutyl) + radical(Propargyl)"""),
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
    label = 'C=[C][CH]C=C(15984)',
    structure = SMILES('[CH2]C=C[C]=C'),
    E0 = (374.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(2.0825,'amu*angstrom^2'), symmetry=1, barrier=(47.8809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07959,'amu*angstrom^2'), symmetry=1, barrier=(47.8139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2251,0.0302874,1.05284e-05,-3.68284e-08,1.72731e-11,45167.7,17.3268], Tmin=(100,'K'), Tmax=(923.512,'K')), NASAPolynomial(coeffs=[10.3879,0.0174193,-5.09535e-06,8.16559e-10,-5.51188e-14,42701,-26.5962], Tmin=(923.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
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
    label = '[CH2][C]=C[CH][CH2](16919)',
    structure = SMILES('[CH2][C]=C[CH][CH2]'),
    E0 = (683.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1114.42],'cm^-1')),
        HinderedRotor(inertia=(0.00740552,'amu*angstrom^2'), symmetry=1, barrier=(6.52457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.14018,'amu*angstrom^2'), symmetry=1, barrier=(95.191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0340784,'amu*angstrom^2'), symmetry=1, barrier=(30.0845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16721,0.0338824,-8.66428e-06,-9.5211e-09,5.10585e-12,82235.7,20.5922], Tmin=(100,'K'), Tmax=(1070.03,'K')), NASAPolynomial(coeffs=[9.10766,0.0208121,-8.38972e-06,1.55219e-09,-1.08387e-13,80013.4,-16.8053], Tmin=(1070.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]([CH2])C=C(17280)',
    structure = SMILES('[CH2][C]=CC([CH2])=C[CH2]'),
    E0 = (609.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,562.923],'cm^-1')),
        HinderedRotor(inertia=(0.284772,'amu*angstrom^2'), symmetry=1, barrier=(63.8987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284826,'amu*angstrom^2'), symmetry=1, barrier=(63.8368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284774,'amu*angstrom^2'), symmetry=1, barrier=(63.8838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283281,'amu*angstrom^2'), symmetry=1, barrier=(63.843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05238,0.0526186,-4.8284e-06,-3.37637e-08,1.77357e-11,73451,26.6575], Tmin=(100,'K'), Tmax=(952.27,'K')), NASAPolynomial(coeffs=[14.6495,0.0251847,-8.36739e-06,1.44442e-09,-1.00284e-13,69515.6,-45.3419], Tmin=(952.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC([CH2])[C]=C(17281)',
    structure = SMILES('[CH2][C]=CC([CH2])[C]=C'),
    E0 = (851.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1670,1700,300,440,367.224,367.423],'cm^-1')),
        HinderedRotor(inertia=(0.0746498,'amu*angstrom^2'), symmetry=1, barrier=(7.14554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746324,'amu*angstrom^2'), symmetry=1, barrier=(7.14599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746136,'amu*angstrom^2'), symmetry=1, barrier=(7.14575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.653593,'amu*angstrom^2'), symmetry=1, barrier=(62.5797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01168,0.0710235,-7.38545e-05,3.90021e-08,-4.86818e-12,102505,27.9886], Tmin=(100,'K'), Tmax=(649.865,'K')), NASAPolynomial(coeffs=[8.45328,0.0353862,-1.50639e-05,2.76471e-09,-1.88633e-13,101323,-6.35725], Tmin=(649.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(851.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])C=C(17282)',
    structure = SMILES('[CH2][C]=[C]C([CH2])C=C'),
    E0 = (851.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1670,1700,300,440,367.224,367.423],'cm^-1')),
        HinderedRotor(inertia=(0.0746498,'amu*angstrom^2'), symmetry=1, barrier=(7.14554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746324,'amu*angstrom^2'), symmetry=1, barrier=(7.14599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746136,'amu*angstrom^2'), symmetry=1, barrier=(7.14575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.653593,'amu*angstrom^2'), symmetry=1, barrier=(62.5797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01168,0.0710235,-7.38545e-05,3.90021e-08,-4.86818e-12,102505,27.9886], Tmin=(100,'K'), Tmax=(649.865,'K')), NASAPolynomial(coeffs=[8.45328,0.0353862,-1.50639e-05,2.76471e-09,-1.88633e-13,101323,-6.35725], Tmin=(649.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(851.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])[CH][C]=C(16337)',
    structure = SMILES('[CH]=CC([CH2])C=[C][CH2]'),
    E0 = (860.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.6502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155481,'amu*angstrom^2'), symmetry=1, barrier=(3.57481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713687,'amu*angstrom^2'), symmetry=1, barrier=(16.4091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86966,'amu*angstrom^2'), symmetry=1, barrier=(65.9791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797695,0.0724935,-7.75867e-05,4.82855e-08,-1.24602e-11,103629,28.6571], Tmin=(100,'K'), Tmax=(931.974,'K')), NASAPolynomial(coeffs=[10.0712,0.0326914,-1.3525e-05,2.45989e-09,-1.67489e-13,101900,-15.429], Tmin=(931.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=C[C](C)C=C(17283)',
    structure = SMILES('[CH2][C]=CC(C)=C[CH2]'),
    E0 = (458.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02644,0.0554289,-1.45321e-05,-1.88232e-08,1.10043e-11,55228.6,26.4698], Tmin=(100,'K'), Tmax=(984.211,'K')), NASAPolynomial(coeffs=[12.5592,0.0308958,-1.11866e-05,1.97101e-09,-1.35256e-13,51876.5,-34.4817], Tmin=(984.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]C([CH2])C=C(17284)',
    structure = SMILES('[CH2]C=[C]C([CH2])C=C'),
    E0 = (613.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,269.17,3233.31],'cm^-1')),
        HinderedRotor(inertia=(0.259691,'amu*angstrom^2'), symmetry=1, barrier=(13.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5476,'amu*angstrom^2'), symmetry=1, barrier=(79.5686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259692,'amu*angstrom^2'), symmetry=1, barrier=(13.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54761,'amu*angstrom^2'), symmetry=1, barrier=(79.5686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882232,0.0689034,-6.09384e-05,3.11838e-08,-6.748e-12,73908.3,27.9463], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[10.0801,0.0350785,-1.42922e-05,2.59388e-09,-1.76871e-13,71907.4,-17.2016], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC(C)[C]=C(17285)',
    structure = SMILES('[CH2][C]=CC(C)[C]=C'),
    E0 = (646.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,317.624,317.644],'cm^-1')),
        HinderedRotor(inertia=(0.129903,'amu*angstrom^2'), symmetry=1, barrier=(9.29444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12969,'amu*angstrom^2'), symmetry=1, barrier=(9.29245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12984,'amu*angstrom^2'), symmetry=1, barrier=(9.29455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367226,'amu*angstrom^2'), symmetry=1, barrier=(26.232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0234,0.0694176,-6.42663e-05,3.5675e-08,-8.6339e-12,77840.3,26.1417], Tmin=(100,'K'), Tmax=(962.377,'K')), NASAPolynomial(coeffs=[8.22964,0.0394661,-1.75831e-05,3.3365e-09,-2.33313e-13,76453.2,-8.34804], Tmin=(962.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(C)C=C(17286)',
    structure = SMILES('[CH2][C]=[C]C(C)C=C'),
    E0 = (646.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,317.624,317.644],'cm^-1')),
        HinderedRotor(inertia=(0.129903,'amu*angstrom^2'), symmetry=1, barrier=(9.29444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12969,'amu*angstrom^2'), symmetry=1, barrier=(9.29245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12984,'amu*angstrom^2'), symmetry=1, barrier=(9.29455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367226,'amu*angstrom^2'), symmetry=1, barrier=(26.232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0234,0.0694176,-6.42663e-05,3.5675e-08,-8.6339e-12,77840.3,26.1417], Tmin=(100,'K'), Tmax=(962.377,'K')), NASAPolynomial(coeffs=[8.22964,0.0394661,-1.75831e-05,3.3365e-09,-2.33313e-13,76453.2,-8.34804], Tmin=(962.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[C]([CH2])C=C(17287)',
    structure = SMILES('[CH2]C=CC([CH2])=C[CH2]'),
    E0 = (371.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0849,0.0480202,1.99814e-05,-6.31139e-08,2.87651e-11,44847.7,26.0673], Tmin=(100,'K'), Tmax=(946.155,'K')), NASAPolynomial(coeffs=[15.8402,0.0256832,-8.08925e-06,1.39535e-09,-9.89541e-14,40263.1,-53.7737], Tmin=(946.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([C]=[C]C)C=C(17288)',
    structure = SMILES('[CH2]C([C]=[C]C)C=C'),
    E0 = (699.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,260.052,260.06],'cm^-1')),
        HinderedRotor(inertia=(0.153535,'amu*angstrom^2'), symmetry=1, barrier=(7.36763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153525,'amu*angstrom^2'), symmetry=1, barrier=(7.36762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(7.36769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153532,'amu*angstrom^2'), symmetry=1, barrier=(7.36767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741746,0.0773159,-9.92364e-05,8.10191e-08,-2.72757e-11,84292.9,28.6415], Tmin=(100,'K'), Tmax=(837.833,'K')), NASAPolynomial(coeffs=[6.48663,0.0408329,-1.7707e-05,3.2455e-09,-2.19537e-13,83648.1,3.83907], Tmin=(837.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC([CH2])[C]=C(17289)',
    structure = SMILES('[CH2]C=CC([CH2])[C]=C'),
    E0 = (613.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,269.17,3233.31],'cm^-1')),
        HinderedRotor(inertia=(0.259691,'amu*angstrom^2'), symmetry=1, barrier=(13.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5476,'amu*angstrom^2'), symmetry=1, barrier=(79.5686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259692,'amu*angstrom^2'), symmetry=1, barrier=(13.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54761,'amu*angstrom^2'), symmetry=1, barrier=(79.5686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882232,0.0689034,-6.09384e-05,3.11838e-08,-6.748e-12,73908.3,27.9463], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[10.0801,0.0350785,-1.42922e-05,2.59388e-09,-1.76871e-13,71907.4,-17.2016], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)[CH][C]=C(16341)',
    structure = SMILES('[CH]=CC(C)C=[C][CH2]'),
    E0 = (655.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,267.237],'cm^-1')),
        HinderedRotor(inertia=(0.21784,'amu*angstrom^2'), symmetry=1, barrier=(11.0458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217679,'amu*angstrom^2'), symmetry=1, barrier=(11.044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217903,'amu*angstrom^2'), symmetry=1, barrier=(11.0456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547482,'amu*angstrom^2'), symmetry=1, barrier=(27.7635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925441,0.0691458,-5.94881e-05,2.87375e-08,-5.89394e-12,78958.8,26.413], Tmin=(100,'K'), Tmax=(1134.19,'K')), NASAPolynomial(coeffs=[10.343,0.0359323,-1.55621e-05,2.91818e-09,-2.0278e-13,76822.6,-20.2071], Tmin=(1134.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](C=C)C=[C]C(17290)',
    structure = SMILES('[CH2]C=C([CH2])C=[C]C'),
    E0 = (491.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.737026,0.0619373,-2.95434e-05,-5.90172e-09,6.89005e-12,59261.1,26.2315], Tmin=(100,'K'), Tmax=(1002.85,'K')), NASAPolynomial(coeffs=[14.1779,0.0288104,-1.06325e-05,1.89409e-09,-1.30651e-13,55535.2,-43.7869], Tmin=(1002.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C=C[CH2](14044)',
    structure = SMILES('[CH]=CC([CH2])C=C[CH2]'),
    E0 = (622.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,243.089],'cm^-1')),
        HinderedRotor(inertia=(0.121276,'amu*angstrom^2'), symmetry=1, barrier=(5.08499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121252,'amu*angstrom^2'), symmetry=1, barrier=(5.08511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.511234,'amu*angstrom^2'), symmetry=1, barrier=(21.4379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.71109,'amu*angstrom^2'), symmetry=1, barrier=(71.747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683213,0.069731,-5.95993e-05,2.8203e-08,-5.51049e-12,75031.5,28.5865], Tmin=(100,'K'), Tmax=(1213.4,'K')), NASAPolynomial(coeffs=[12.3135,0.0313918,-1.22049e-05,2.16374e-09,-1.45602e-13,72209,-29.7725], Tmin=(1213.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([C]=C)C=[C]C(17291)',
    structure = SMILES('[CH2]C([C]=C)C=[C]C'),
    E0 = (699.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,260.052,260.06],'cm^-1')),
        HinderedRotor(inertia=(0.153535,'amu*angstrom^2'), symmetry=1, barrier=(7.36763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153525,'amu*angstrom^2'), symmetry=1, barrier=(7.36762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153521,'amu*angstrom^2'), symmetry=1, barrier=(7.36769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153532,'amu*angstrom^2'), symmetry=1, barrier=(7.36767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741746,0.0773159,-9.92364e-05,8.10191e-08,-2.72757e-11,84292.9,28.6415], Tmin=(100,'K'), Tmax=(837.833,'K')), NASAPolynomial(coeffs=[6.48663,0.0408329,-1.7707e-05,3.2455e-09,-2.19537e-13,83648.1,3.83907], Tmin=(837.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C=[C]C(17292)',
    structure = SMILES('[CH]=CC([CH2])C=[C]C'),
    E0 = (709.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.315048,'amu*angstrom^2'), symmetry=1, barrier=(7.24357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314779,'amu*angstrom^2'), symmetry=1, barrier=(7.23739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.315372,'amu*angstrom^2'), symmetry=1, barrier=(7.25101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577576,'amu*angstrom^2'), symmetry=1, barrier=(13.2796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13234,0.0699565,-6.18508e-05,1.60867e-08,1.04481e-11,85390.8,27.2365], Tmin=(100,'K'), Tmax=(588.469,'K')), NASAPolynomial(coeffs=[7.61611,0.0390746,-1.67549e-05,3.08771e-09,-2.11131e-13,84399.3,-2.54659], Tmin=(588.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C([C]=C)C=C(17293)',
    structure = SMILES('[CH2][CH]C([C]=C)C=C'),
    E0 = (669.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,346.942,356.216,1313.27],'cm^-1')),
        HinderedRotor(inertia=(0.00137082,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853597,'amu*angstrom^2'), symmetry=1, barrier=(7.3913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858293,'amu*angstrom^2'), symmetry=1, barrier=(7.37088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0838535,'amu*angstrom^2'), symmetry=1, barrier=(7.3618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17927,0.0672513,-6.73655e-05,4.54645e-08,-1.39939e-11,80562.6,29.6122], Tmin=(100,'K'), Tmax=(757.817,'K')), NASAPolynomial(coeffs=[5.70529,0.0433624,-2.00822e-05,3.86988e-09,-2.72473e-13,79876.6,9.0318], Tmin=(757.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1CC1C=C(17294)',
    structure = SMILES('C=[C]C1CC1C=C'),
    E0 = (424.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49688,0.0363406,4.54545e-05,-8.57935e-08,3.56126e-11,51152.3,25.1133], Tmin=(100,'K'), Tmax=(961.991,'K')), NASAPolynomial(coeffs=[15.5667,0.0242667,-8.11443e-06,1.50061e-09,-1.11408e-13,46296.9,-53.3863], Tmin=(961.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(=C)C=C(15829)',
    structure = SMILES('C=[C]CC(=C)C=C'),
    E0 = (376.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,1685,370,221.256,221.66,222.606],'cm^-1')),
        HinderedRotor(inertia=(0.490568,'amu*angstrom^2'), symmetry=1, barrier=(17.0028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486786,'amu*angstrom^2'), symmetry=1, barrier=(17.0024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488039,'amu*angstrom^2'), symmetry=1, barrier=(17.0028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712864,0.0615411,-2.86001e-05,-8.26196e-09,7.99054e-12,45446.2,26.526], Tmin=(100,'K'), Tmax=(1004.13,'K')), NASAPolynomial(coeffs=[15.3147,0.0261308,-9.69752e-06,1.75775e-09,-1.2322e-13,41366.5,-49.6923], Tmin=(1004.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C(C)C=C(17295)',
    structure = SMILES('[CH2]C=C(C)C=C=C'),
    E0 = (278.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816619,0.0593926,-2.20705e-05,-1.3735e-08,9.66525e-12,33672.7,24.2453], Tmin=(100,'K'), Tmax=(994.941,'K')), NASAPolynomial(coeffs=[14.2369,0.0285748,-1.04895e-05,1.87716e-09,-1.30364e-13,29857,-46.1865], Tmin=(994.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C1[CH]CC1[C]=C(17296)',
    structure = SMILES('[CH2]C1[CH]CC1[C]=C'),
    E0 = (684.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65896,0.0376596,2.99247e-05,-6.10887e-08,2.4948e-11,82370.5,28.1663], Tmin=(100,'K'), Tmax=(978.419,'K')), NASAPolynomial(coeffs=[11.5095,0.0309596,-1.12713e-05,2.04959e-09,-1.45219e-13,78836,-27.3536], Tmin=(978.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1C([CH2])C1[C]=C(17297)',
    structure = SMILES('[CH2]C1C([CH2])C1[C]=C'),
    E0 = (696.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29124,0.0445947,2.19574e-05,-6.50503e-08,3.02682e-11,83903.3,27.9343], Tmin=(100,'K'), Tmax=(920.671,'K')), NASAPolynomial(coeffs=[15.185,0.022992,-5.99779e-06,9.20805e-10,-6.2817e-14,79702.2,-46.8683], Tmin=(920.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C([CH2])C=C(17298)',
    structure = SMILES('C#C[CH]C([CH2])C=C'),
    E0 = (549.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,391.815,2078.83],'cm^-1')),
        HinderedRotor(inertia=(0.664652,'amu*angstrom^2'), symmetry=1, barrier=(72.34,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662956,'amu*angstrom^2'), symmetry=1, barrier=(72.3229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0898732,'amu*angstrom^2'), symmetry=1, barrier=(9.81269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662721,'amu*angstrom^2'), symmetry=1, barrier=(72.3217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420631,0.0672264,-5.98011e-05,2.93022e-08,-5.63479e-12,66209.5,29.3489], Tmin=(100,'K'), Tmax=(1415.18,'K')), NASAPolynomial(coeffs=[14.0211,0.0225391,-5.81531e-06,7.5183e-10,-4.02747e-14,62985.5,-38.7785], Tmin=(1415.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C=CC=C(17299)',
    structure = SMILES('[CH2]C=CC=C=C'),
    E0 = (316.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.57103,'amu*angstrom^2'), symmetry=1, barrier=(36.1211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57001,'amu*angstrom^2'), symmetry=1, barrier=(36.0977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67668,0.0391992,9.13492e-06,-4.16224e-08,1.94021e-11,38164.2,20.5575], Tmin=(100,'K'), Tmax=(956.854,'K')), NASAPolynomial(coeffs=[13.1804,0.0202032,-6.69471e-06,1.18319e-09,-8.42505e-14,34630.8,-41.3935], Tmin=(956.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C](C=C)C[C]=C(15832)',
    structure = SMILES('[CH2]C=C([CH2])C[C]=C'),
    E0 = (554.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,809.622,810.159],'cm^-1')),
        HinderedRotor(inertia=(0.0792804,'amu*angstrom^2'), symmetry=1, barrier=(1.82281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896671,'amu*angstrom^2'), symmetry=1, barrier=(20.6162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896774,'amu*angstrom^2'), symmetry=1, barrier=(20.6186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51147,'amu*angstrom^2'), symmetry=1, barrier=(34.7516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650314,0.0639325,-3.69869e-05,3.39338e-09,2.96099e-12,66821,28.447], Tmin=(100,'K'), Tmax=(1072.01,'K')), NASAPolynomial(coeffs=[14.5779,0.0287818,-1.13346e-05,2.0745e-09,-1.44151e-13,62868.5,-44.2214], Tmin=(1072.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C([CH2])C=C(14210)',
    structure = SMILES('[CH]C=CC([CH2])C=C'),
    E0 = (594.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854398,0.067509,-4.72676e-05,1.86115e-08,-3.173e-12,71666.9,28.3946], Tmin=(100,'K'), Tmax=(1320.6,'K')), NASAPolynomial(coeffs=[9.91908,0.0400527,-1.60815e-05,2.86811e-09,-1.92645e-13,69272.7,-17.858], Tmin=(1320.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(594.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([C]=C)C[C]=C(15834)',
    structure = SMILES('[CH2]C([C]=C)C[C]=C'),
    E0 = (712.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0305923,'amu*angstrom^2'), symmetry=1, barrier=(12.5243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110226,'amu*angstrom^2'), symmetry=1, barrier=(12.5238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00204458,'amu*angstrom^2'), symmetry=1, barrier=(12.5236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192736,'amu*angstrom^2'), symmetry=1, barrier=(78.889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951852,0.0653319,-5.31855e-05,2.44992e-08,-4.75174e-12,85830.8,30.71], Tmin=(100,'K'), Tmax=(1205.64,'K')), NASAPolynomial(coeffs=[10.656,0.0331362,-1.31292e-05,2.34976e-09,-1.58869e-13,83490.8,-17.9215], Tmin=(1205.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(712.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC([CH2])C=C(15835)',
    structure = SMILES('[CH]=[C]CC([CH2])C=C'),
    E0 = (721.968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,1089.63],'cm^-1')),
        HinderedRotor(inertia=(0.695087,'amu*angstrom^2'), symmetry=1, barrier=(15.9814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00246626,'amu*angstrom^2'), symmetry=1, barrier=(15.9831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566213,'amu*angstrom^2'), symmetry=1, barrier=(13.0184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.56089,'amu*angstrom^2'), symmetry=1, barrier=(81.8719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693257,0.0667862,-5.3715e-05,2.35438e-08,-4.23194e-12,86956.7,31.5694], Tmin=(100,'K'), Tmax=(1318.34,'K')), NASAPolynomial(coeffs=[13.0986,0.0291467,-1.08889e-05,1.88721e-09,-1.25141e-13,83685.8,-31.7079], Tmin=(1318.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(721.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C[C]=C(15590)',
    structure = SMILES('[CH]=CC([CH2])C[C]=C'),
    E0 = (721.968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,1089.63],'cm^-1')),
        HinderedRotor(inertia=(0.695087,'amu*angstrom^2'), symmetry=1, barrier=(15.9814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00246626,'amu*angstrom^2'), symmetry=1, barrier=(15.9831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566213,'amu*angstrom^2'), symmetry=1, barrier=(13.0184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.56089,'amu*angstrom^2'), symmetry=1, barrier=(81.8719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3375.81,'J/mol'), sigma=(5.99702,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.29 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693257,0.0667862,-5.3715e-05,2.35438e-08,-4.23194e-12,86956.7,31.5694], Tmin=(100,'K'), Tmax=(1318.34,'K')), NASAPolynomial(coeffs=[13.0986,0.0291467,-1.08889e-05,1.88721e-09,-1.25141e-13,83685.8,-31.7079], Tmin=(1318.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(721.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(C)C=C(17300)',
    structure = SMILES('[CH][C]=CC(C)C=C'),
    E0 = (627.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12196,0.0666526,-4.63336e-05,1.8244e-08,-3.23535e-12,75593.1,26.1284], Tmin=(100,'K'), Tmax=(1223.01,'K')), NASAPolynomial(coeffs=[7.79672,0.0448219,-1.95586e-05,3.64881e-09,-2.51875e-13,73960.5,-7.41714], Tmin=(1223.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    E0 = (613.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (773.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (773.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (933.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1068.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1044.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (621.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (676.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (846.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (992.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (943.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (737.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (651.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (701.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (650.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (667.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (733.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (654.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (742.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (670.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (732.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (668.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (970.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (889.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (821.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1063.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1063.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1072.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (712.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (817.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (798.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (798.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (806.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (861.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (714.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (699.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (744.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (655.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (842.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (841.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (832.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (618.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (636.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (636.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (753.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (696.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (776.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (737.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (771.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (808.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (854.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (867.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (766.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (660.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['C3H3(5450)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2][C]=C[CH]CC=C(17265)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C=CC[CH][C]=C(15872)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH2][C]=CC=C[CH2](17118)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]C(C=C)C=[C][CH2](17266)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH][C]=CC([CH2])C=C(17267)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1=CC(C=C)C1(17268)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C=CC(=C)C=C(17269)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=C[C]([CH2])C[CH2](17270)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=[C]C([CH2])C[CH2](17271)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=[C]C([CH2])[CH]C(17272)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['C=[C][CH]C1[CH]CC1(16131)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1=CC([CH2])[CH]C1(17273)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.80657e+08,'s^-1'), n=0.835, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_pri_2H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1[CH]CC[C]=C1(17274)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2][C]=CC1CC1[CH2](17275)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1=CC([CH2])C1[CH2](17276)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_2H_pri;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1C=[C]CC1[CH2](17277)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][C]=CC(=C)C=C(17278)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-TwoDe_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2]C#CC([CH2])C=C(17279)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H3(60)', 'C=[C][CH]C=C(15984)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00802013,'m^3/(mol*s)'), n=2.41, Ea=(8.77571,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C][CH2](16918)', 'butadiene13(1350)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C3H3(5450)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.08e-06,'m^3/(mol*s)'), n=3.05, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CdsJ=Cdd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C2H3(60)', '[CH2][C]=C[CH][CH2](16919)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.4834e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C][CH2](16918)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][C]=C[C]([CH2])C=C(17280)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2][C]=CC([CH2])[C]=C(17281)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][C]=[C]C([CH2])C=C(17282)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2][C]=C[C](C)C=C(17283)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C=[C]C([CH2])C=C(17284)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=CC(C)[C]=C(17285)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]=[C]C(C)C=C(17286)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C=C[C]([CH2])C=C(17287)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([C]=[C]C)C=C(17288)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=CC([CH2])[C]=C(17289)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSD;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=CC(C)[CH][C]=C(16341)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2][C](C=C)C=[C]C(17290)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC([CH2])C=C[CH2](14044)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([C]=C)C=[C]C(17291)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=CC([CH2])C=[C]C(17292)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_4;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([C]=C)C=C(17293)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['C=[C]C1CC1C=C(17294)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['C=[C]CC(=C)C=C(15829)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['C=[C]C=C(C)C=C(17295)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1[CH]CC1[C]=C(17296)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.61918e+08,'s^-1'), n=0.930343, Ea=(139.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2]C1C([CH2])C1[C]=C(17297)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(9.65211e+08,'s^-1'), n=0.77, Ea=(83.1077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 79.8 to 83.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction47',
    reactants = ['H(3)', 'C#C[CH]C([CH2])C=C(17298)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CH2(T)(20)', 'C=[C]C=CC=C(17299)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH2][C](C=C)C[C]=C(15832)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([CH][C]=C)C=C(15833)'],
    products = ['[CH]=C[CH]C([CH2])C=C(14210)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C([C]=C)C[C]=C(15834)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=[C]CC([CH2])C=C(15835)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=CC([CH2])C[C]=C(15590)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]=[C][CH]C(C)C=C(17300)'],
    products = ['[CH2]C([CH][C]=C)C=C(15833)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3965',
    isomers = [
        '[CH2]C([CH][C]=C)C=C(15833)',
    ],
    reactants = [
        ('C3H3(5450)', 'butadiene13(1350)'),
        ('C3H3(5450)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3965',
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


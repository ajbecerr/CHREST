species(
    label = '[CH]=C(C)C[CH]C(633)',
    structure = SMILES('[CH]=C(C)C[CH]C'),
    E0 = (362.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.123464,'amu*angstrom^2'), symmetry=1, barrier=(2.83868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123485,'amu*angstrom^2'), symmetry=1, barrier=(2.83916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12348,'amu*angstrom^2'), symmetry=1, barrier=(2.83904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123475,'amu*angstrom^2'), symmetry=1, barrier=(2.83894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69405,0.0540655,-3.34381e-05,1.15001e-08,-1.80015e-12,43640.7,24.4885], Tmin=(100,'K'), Tmax=(1323.25,'K')), NASAPolynomial(coeffs=[6.56268,0.0393483,-1.67551e-05,3.095e-09,-2.12187e-13,42352.2,-0.363434], Tmin=(1323.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = 'C#CC(317)',
    structure = SMILES('C#CC'),
    E0 = (172.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.58471,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2748.36,'J/mol'), sigma=(4.8439,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=429.29 K, Pc=54.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30775,0.010871,1.58532e-05,-2.85512e-08,1.24304e-11,20724.7,7.26856], Tmin=(100,'K'), Tmax=(883.298,'K')), NASAPolynomial(coeffs=[6.14675,0.00980871,-2.37135e-06,3.20148e-10,-1.94935e-14,19763.1,-8.6801], Tmin=(883.298,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C=CC(38)',
    structure = SMILES('C=CC'),
    E0 = (6.12372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.597443,'amu*angstrom^2'), symmetry=1, barrier=(13.7364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30973,0.00827538,3.377e-05,-4.39288e-08,1.58764e-11,767.477,9.64363], Tmin=(100,'K'), Tmax=(988.015,'K')), NASAPolynomial(coeffs=[5.4122,0.0172863,-6.51344e-06,1.20319e-09,-8.55894e-14,-503.241,-4.80239], Tmin=(988.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.12372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CH2(S)(2)',
    structure = SMILES('[CH2]'),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144068,5.45069e-06,-3.58002e-09,7.56192e-13,50400.6,-0.411765], Tmin=(100,'K'), Tmax=(1442.36,'K')), NASAPolynomial(coeffs=[2.62647,0.00394763,-1.49924e-06,2.54539e-10,-1.62956e-14,50691.8,6.78379], Tmin=(1442.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=CC[CH]C(328)',
    structure = SMILES('[CH]=CC[CH]C'),
    E0 = (401.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0045432,'amu*angstrom^2'), symmetry=1, barrier=(5.15925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225032,'amu*angstrom^2'), symmetry=1, barrier=(5.17392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227418,'amu*angstrom^2'), symmetry=1, barrier=(5.22878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3052.11,'J/mol'), sigma=(5.53315,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.73 K, Pc=40.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05382,0.0394038,-2.00486e-05,4.87108e-09,-4.71721e-13,48331,22.4667], Tmin=(100,'K'), Tmax=(2304.47,'K')), NASAPolynomial(coeffs=[13.268,0.0199387,-7.37853e-06,1.2057e-09,-7.40813e-14,43162.5,-40.9972], Tmin=(2304.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C)C([CH2])C(635)',
    structure = SMILES('[CH]=C(C)C([CH2])C'),
    E0 = (369.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0930987,'amu*angstrom^2'), symmetry=1, barrier=(5.22619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0930924,'amu*angstrom^2'), symmetry=1, barrier=(5.22626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0931002,'amu*angstrom^2'), symmetry=1, barrier=(5.22632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385515,'amu*angstrom^2'), symmetry=1, barrier=(21.6434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.85,'J/mol'), sigma=(5.90911,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.46 K, Pc=35.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01976,0.0586178,-3.51251e-05,6.72555e-09,1.23172e-12,44497.7,25.7158], Tmin=(100,'K'), Tmax=(1046.42,'K')), NASAPolynomial(coeffs=[11.4005,0.0304178,-1.11593e-05,1.94212e-09,-1.30434e-13,41696.6,-27.8401], Tmin=(1046.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(16)',
    structure = SMILES('[CH]C'),
    E0 = (351.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,431.536,1804.51],'cm^-1')),
        HinderedRotor(inertia=(0.000906354,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73284,-0.0002446,3.59192e-05,-4.4428e-08,1.65883e-11,42287.5,7.07805], Tmin=(100,'K'), Tmax=(940.484,'K')), NASAPolynomial(coeffs=[5.42973,0.00816763,-2.42525e-06,4.22631e-10,-3.09408e-14,41277.1,-4.67916], Tmin=(940.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C(557)',
    structure = SMILES('[CH]C(=C)C'),
    E0 = (337.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,296.624,296.937,297.046],'cm^-1')),
        HinderedRotor(inertia=(0.816008,'amu*angstrom^2'), symmetry=1, barrier=(50.9344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816462,'amu*angstrom^2'), symmetry=1, barrier=(50.9278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45139,0.0283382,3.33456e-06,-1.74484e-08,6.92052e-12,40682.9,14.5918], Tmin=(100,'K'), Tmax=(1068.99,'K')), NASAPolynomial(coeffs=[6.47597,0.0256533,-1.0262e-05,1.8599e-09,-1.2757e-13,39115.4,-8.40003], Tmin=(1068.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH3](12)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.567,1408.02,1408.6,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48743e-06,-3.32749e-09,8.49963e-13,16285.6,0.35174], Tmin=(100,'K'), Tmax=(1337.62,'K')), NASAPolynomial(coeffs=[3.54145,0.00476788,-1.82149e-06,3.28878e-10,-2.22547e-14,16224,1.66039], Tmin=(1337.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH]CC(=[CH])C(558)',
    structure = SMILES('[CH]CC(=[CH])C'),
    E0 = (639.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,180,1318.34,2406.82],'cm^-1')),
        HinderedRotor(inertia=(0.0662821,'amu*angstrom^2'), symmetry=1, barrier=(16.9421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220732,'amu*angstrom^2'), symmetry=1, barrier=(5.07506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.735222,'amu*angstrom^2'), symmetry=1, barrier=(16.9042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67621,0.0468505,-3.55013e-05,1.42174e-08,-2.33569e-12,77029.4,20.5857], Tmin=(100,'K'), Tmax=(1424.54,'K')), NASAPolynomial(coeffs=[10.8874,0.020986,-8.26646e-06,1.47172e-09,-9.88669e-14,74405.1,-27.1122], Tmin=(1424.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = 'H(9)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C(C)C[C]C(704)',
    structure = SMILES('[CH]=C(C)C[C]C'),
    E0 = (615.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00882422,'amu*angstrom^2'), symmetry=1, barrier=(9.02301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772366,'amu*angstrom^2'), symmetry=1, barrier=(17.7582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76801,'amu*angstrom^2'), symmetry=1, barrier=(17.6581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.572924,'amu*angstrom^2'), symmetry=1, barrier=(13.1726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10449,0.0608123,-4.67305e-05,1.92025e-08,-3.27805e-12,74188.6,23.7804], Tmin=(100,'K'), Tmax=(1359.05,'K')), NASAPolynomial(coeffs=[11.776,0.0294028,-1.20627e-05,2.19626e-09,-1.4964e-13,71288.1,-30.977], Tmin=(1359.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=C(C)C[CH]C(705)',
    structure = SMILES('[C]=C(C)C[CH]C'),
    E0 = (673.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,350,440,435,1725,1647.01,1647.04],'cm^-1')),
        HinderedRotor(inertia=(0.14988,'amu*angstrom^2'), symmetry=1, barrier=(8.43603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150328,'amu*angstrom^2'), symmetry=1, barrier=(8.43512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149989,'amu*angstrom^2'), symmetry=1, barrier=(8.43517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149917,'amu*angstrom^2'), symmetry=1, barrier=(8.43555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48864,0.0596827,-6.35006e-05,5.04001e-08,-1.8097e-11,81052.4,25.21], Tmin=(100,'K'), Tmax=(765.778,'K')), NASAPolynomial(coeffs=[4.12546,0.0413354,-1.86024e-05,3.5129e-09,-2.43483e-13,80782.7,14.0682], Tmin=(765.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'CC1=CC(C)C1(695)',
    structure = SMILES('CC1=CC(C)C1'),
    E0 = (73.6604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87849,0.0290293,5.6384e-05,-9.03629e-08,3.56451e-11,8951.51,18.9778], Tmin=(100,'K'), Tmax=(970.034,'K')), NASAPolynomial(coeffs=[12.7505,0.0280649,-9.95804e-06,1.85058e-09,-1.35256e-13,4778.38,-43.7809], Tmin=(970.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.6604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C(C)C=CC(677)',
    structure = SMILES('C=C(C)C=CC'),
    E0 = (20.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14009,0.0531094,-1.57896e-05,-1.45952e-08,8.71461e-12,2553.6,20.6675], Tmin=(100,'K'), Tmax=(1023.2,'K')), NASAPolynomial(coeffs=[12.4481,0.0294413,-1.12018e-05,2.0334e-09,-1.41508e-13,-835.616,-39.4005], Tmin=(1023.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(=C)C(629)',
    structure = SMILES('C=CCC(=C)C'),
    E0 = (49.6948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35878,0.048418,-7.65502e-06,-1.88636e-08,9.12828e-12,6080.25,23.0633], Tmin=(100,'K'), Tmax=(1056.43,'K')), NASAPolynomial(coeffs=[11.0985,0.0317129,-1.25786e-05,2.31868e-09,-1.61897e-13,2896.72,-29.7874], Tmin=(1056.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.6948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=C(C)C=CC(706)',
    structure = SMILES('[CH]=C(C)C=CC'),
    E0 = (267.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.728206,'amu*angstrom^2'), symmetry=1, barrier=(16.7429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727708,'amu*angstrom^2'), symmetry=1, barrier=(16.7314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728418,'amu*angstrom^2'), symmetry=1, barrier=(16.7478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04677,0.0568551,-3.31855e-05,3.67197e-09,2.45803e-12,32274.2,21.406], Tmin=(100,'K'), Tmax=(1054.82,'K')), NASAPolynomial(coeffs=[12.8292,0.0263754,-1.00361e-05,1.80397e-09,-1.24125e-13,28998.5,-39.8112], Tmin=(1054.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C)CC=C(707)',
    structure = SMILES('[CH]=C(C)CC=C'),
    E0 = (296.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,303.373],'cm^-1')),
        HinderedRotor(inertia=(0.185512,'amu*angstrom^2'), symmetry=1, barrier=(12.0959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185942,'amu*angstrom^2'), symmetry=1, barrier=(12.0968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185001,'amu*angstrom^2'), symmetry=1, barrier=(12.0902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25732,0.0522413,-2.52311e-05,-4.99264e-10,2.89078e-12,35801.2,23.8322], Tmin=(100,'K'), Tmax=(1111.52,'K')), NASAPolynomial(coeffs=[11.6326,0.0284034,-1.12789e-05,2.05875e-09,-1.42054e-13,32660.9,-31.0706], Tmin=(1111.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(318)',
    structure = SMILES('[CH]=[C]C'),
    E0 = (491.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.144654,'amu*angstrom^2'), symmetry=1, barrier=(3.32589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19643,0.0163496,-6.94933e-06,1.03466e-09,1.39733e-14,59091.1,10.936], Tmin=(100,'K'), Tmax=(1736.56,'K')), NASAPolynomial(coeffs=[6.49591,0.0109619,-4.20665e-06,7.15392e-10,-4.5682e-14,57611.6,-7.7636], Tmin=(1736.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C(39)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (279.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00418542,'amu*angstrom^2'), symmetry=1, barrier=(6.91843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00418543,'amu*angstrom^2'), symmetry=1, barrier=(6.91842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2737.71,'J/mol'), sigma=(5.09232,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=427.62 K, Pc=47.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25505,0.0137285,1.00537e-05,-1.43788e-08,4.38751e-12,33590.4,14.1736], Tmin=(100,'K'), Tmax=(1201.86,'K')), NASAPolynomial(coeffs=[3.74312,0.0203097,-8.40105e-06,1.5386e-09,-1.05137e-13,32880.4,9.26377], Tmin=(1201.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJC) + radical(RCCJ)"""),
)

species(
    label = 'C#CC[CH]C(469)',
    structure = SMILES('C#CC[CH]C'),
    E0 = (320.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.099383,'amu*angstrom^2'), symmetry=1, barrier=(2.28501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0108951,'amu*angstrom^2'), symmetry=1, barrier=(77.8463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.38701,'amu*angstrom^2'), symmetry=1, barrier=(77.8741,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23907,0.0380811,-2.41037e-05,9.4766e-09,-1.71379e-12,38580.8,19.2199], Tmin=(100,'K'), Tmax=(1210.42,'K')), NASAPolynomial(coeffs=[5.6734,0.0267321,-1.00397e-05,1.73069e-09,-1.13976e-13,37749.3,1.99533], Tmin=(1210.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJC)"""),
)

species(
    label = '[CH]=C(C)[CH][CH]C(708)',
    structure = SMILES('[CH]C(C)=C[CH]C'),
    E0 = (420.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,514.876,514.879,514.882,514.882],'cm^-1')),
        HinderedRotor(inertia=(0.279691,'amu*angstrom^2'), symmetry=1, barrier=(52.6162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27969,'amu*angstrom^2'), symmetry=1, barrier=(52.6161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279693,'amu*angstrom^2'), symmetry=1, barrier=(52.6161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27969,'amu*angstrom^2'), symmetry=1, barrier=(52.6161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27461,0.0514161,-1.3273e-05,-1.04131e-08,5.34722e-12,50642.9,22.4326], Tmin=(100,'K'), Tmax=(1140.16,'K')), NASAPolynomial(coeffs=[9.81306,0.0371559,-1.51608e-05,2.76414e-09,-1.89414e-13,47675.7,-24.3542], Tmin=(1140.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[CH]C(471)',
    structure = SMILES('[CH]=[C]C[CH]C'),
    E0 = (639.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,232.443],'cm^-1')),
        HinderedRotor(inertia=(0.00275455,'amu*angstrom^2'), symmetry=1, barrier=(5.88367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153223,'amu*angstrom^2'), symmetry=1, barrier=(5.88522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153684,'amu*angstrom^2'), symmetry=1, barrier=(5.88253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91189,0.0302195,2.7808e-05,-1.26114e-07,1.15622e-10,76896.3,20.0276], Tmin=(100,'K'), Tmax=(409.644,'K')), NASAPolynomial(coeffs=[3.8483,0.0313221,-1.3748e-05,2.57465e-09,-1.78009e-13,76733.6,15.2963], Tmin=(409.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C[CH]C(676)',
    structure = SMILES('[CH]C(=C)C[CH]C'),
    E0 = (485.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4528,0.0539888,-2.77267e-05,6.78589e-09,-6.69264e-13,58517.8,26.0707], Tmin=(100,'K'), Tmax=(2198.3,'K')), NASAPolynomial(coeffs=[14.1168,0.0309455,-1.20032e-05,2.01752e-09,-1.26984e-13,52949.9,-45.0012], Tmin=(2198.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C)C[CH][CH2](709)',
    structure = SMILES('[CH]=C(C)C[CH][CH2]'),
    E0 = (567.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,1576.52],'cm^-1')),
        HinderedRotor(inertia=(0.156177,'amu*angstrom^2'), symmetry=1, barrier=(3.59083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156021,'amu*angstrom^2'), symmetry=1, barrier=(3.58723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156128,'amu*angstrom^2'), symmetry=1, barrier=(3.58968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155958,'amu*angstrom^2'), symmetry=1, barrier=(3.58577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63479,0.0556692,-4.58814e-05,2.52312e-08,-6.53747e-12,68328.3,26.5078], Tmin=(100,'K'), Tmax=(871.541,'K')), NASAPolynomial(coeffs=[5.27714,0.0389523,-1.71098e-05,3.22283e-09,-2.24342e-13,67693.4,9.43645], Tmin=(871.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C)[CH]CC(685)',
    structure = SMILES('[CH]C(C)=CCC'),
    E0 = (279.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,422.345,422.434,422.496,422.552],'cm^-1')),
        HinderedRotor(inertia=(0.418745,'amu*angstrom^2'), symmetry=1, barrier=(53.0316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418533,'amu*angstrom^2'), symmetry=1, barrier=(53.0312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418892,'amu*angstrom^2'), symmetry=1, barrier=(53.0314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418803,'amu*angstrom^2'), symmetry=1, barrier=(53.0322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0272,0.0571387,-2.35742e-05,-7.50338e-10,2.11832e-12,33679.6,24.4326], Tmin=(100,'K'), Tmax=(1238.1,'K')), NASAPolynomial(coeffs=[10.5796,0.0383695,-1.54849e-05,2.78249e-09,-1.87918e-13,30387.4,-27.4352], Tmin=(1238.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C)CC[CH2](710)',
    structure = SMILES('[CH]=C(C)CC[CH2]'),
    E0 = (372.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,767.225],'cm^-1')),
        HinderedRotor(inertia=(0.105804,'amu*angstrom^2'), symmetry=1, barrier=(2.43265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10697,'amu*angstrom^2'), symmetry=1, barrier=(2.45945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.473555,'amu*angstrom^2'), symmetry=1, barrier=(10.888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.473181,'amu*angstrom^2'), symmetry=1, barrier=(10.8794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13585,0.0591304,-3.96319e-05,1.3872e-08,-2.0213e-12,44965.7,25.5512], Tmin=(100,'K'), Tmax=(1553.88,'K')), NASAPolynomial(coeffs=[12.0424,0.0310546,-1.25296e-05,2.2442e-09,-1.50521e-13,41576.3,-31.8736], Tmin=(1553.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(C)[CH][CH]C(683)',
    structure = SMILES('[CH2]C(C)=C[CH]C'),
    E0 = (201.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3144,0.0486071,-4.15497e-06,-2.37794e-08,1.10435e-11,24281.5,21.4565], Tmin=(100,'K'), Tmax=(1043.43,'K')), NASAPolynomial(coeffs=[11.482,0.0321074,-1.27487e-05,2.35694e-09,-1.65195e-13,20936.1,-33.8919], Tmin=(1043.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C[CH]C(436)',
    structure = SMILES('[CH2]C(=C)C[CH]C'),
    E0 = (266.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,397.682,2062.23],'cm^-1')),
        HinderedRotor(inertia=(0.0748507,'amu*angstrom^2'), symmetry=1, barrier=(8.39351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00982772,'amu*angstrom^2'), symmetry=1, barrier=(29.6775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748346,'amu*angstrom^2'), symmetry=1, barrier=(8.39413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854387,'amu*angstrom^2'), symmetry=1, barrier=(96.1175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3224.65,'J/mol'), sigma=(5.87559,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.68 K, Pc=36.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07974,0.0556446,-3.26018e-05,9.62163e-09,-1.15636e-12,32175.6,25.9137], Tmin=(100,'K'), Tmax=(1886.98,'K')), NASAPolynomial(coeffs=[14.3542,0.0275054,-1.02333e-05,1.71883e-09,-1.09339e-13,27165.9,-46.557], Tmin=(1886.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([CH2])CCC(684)',
    structure = SMILES('[CH]C(=C)CCC'),
    E0 = (291.327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01177,0.0566454,-1.84653e-05,-8.19645e-09,5.11569e-12,35153.7,24.9197], Tmin=(100,'K'), Tmax=(1122.26,'K')), NASAPolynomial(coeffs=[10.7458,0.0381242,-1.53273e-05,2.78103e-09,-1.90369e-13,31950.4,-27.7019], Tmin=(1122.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][CH]CC(=C)C(422)',
    structure = SMILES('[CH2][CH]CC(=C)C'),
    E0 = (320.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180,1599.33],'cm^-1')),
        HinderedRotor(inertia=(0.0892705,'amu*angstrom^2'), symmetry=1, barrier=(2.0525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889334,'amu*angstrom^2'), symmetry=1, barrier=(2.04475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0891168,'amu*angstrom^2'), symmetry=1, barrier=(2.04897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889638,'amu*angstrom^2'), symmetry=1, barrier=(2.04545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82264,0.0512649,-2.77873e-05,7.49967e-09,-8.56116e-13,38601.8,25.3925], Tmin=(100,'K'), Tmax=(1779.89,'K')), NASAPolynomial(coeffs=[8.6284,0.0359703,-1.48979e-05,2.67195e-09,-1.78032e-13,36179.1,-11.3654], Tmin=(1779.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
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
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (362.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (820.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (528.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (746.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (831.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (827.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (884.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (370.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (425.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (387.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (487.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (511.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (509.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (480.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (492.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (770.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (631.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (774.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (697.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (779.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (469.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (524.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (554.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (491.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (377.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (395.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['C#CC(317)', 'C=CC(38)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(2)', '[CH]=CC[CH]C(328)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(71881.9,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C)C([CH2])C(635)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(16)', '[CH]=C([CH2])C(557)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH3](12)', '[CH]CC(=[CH])C(558)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(9)', '[CH]=C(C)C[C]C(704)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(9)', '[C]=C(C)C[CH]C(705)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['CC1=CC(C)C1(695)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['C=C(C)C=CC(677)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['C=CCC(=C)C(629)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(9)', '[CH]=C(C)C=CC(706)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(105.529,'m^3/(mol*s)'), n=1.6629, Ea=(8.08712,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(9)', '[CH]=C(C)CC=C(707)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]C(318)', 'C=CC(38)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(12.4666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC(317)', '[CH2][CH]C(39)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0942128,'m^3/(mol*s)'), n=2.31088, Ea=(29.6884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cs_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH3](12)', 'C#CC[CH]C(469)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(138000,'cm^3/(mol*s)'), n=2.41, Ea=(36.8192,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2325 used for Ct-Cs_Ct-H;CsJ-HHH
Exact match found for rate rule [Ct-Cs_Ct-H;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]C(318)', '[CH2][CH]C(39)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(9)', '[CH]=C(C)[CH][CH]C(708)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH3](12)', '[CH]=[C]C[CH]C(471)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.48129e+07,'m^3/(mol*s)'), n=-0.157514, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-3C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-3C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-3C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(9)', '[CH]C(=C)C[CH]C(676)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(9)', '[CH]=C(C)C[CH][CH2](709)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C)[CH]CC(685)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.9e+10,'s^-1'), n=0.75, Ea=(190.79,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 159 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C)CC[CH2](710)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['C=C(C)[CH][CH]C(683)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['[CH2]C(=C)C[CH]C(436)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C([CH2])CCC(684)'],
    products = ['[CH]=C(C)C[CH]C(633)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C)C[CH]C(633)'],
    products = ['[CH2][CH]CC(=C)C(422)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #163',
    isomers = [
        '[CH]=C(C)C[CH]C(633)',
    ],
    reactants = [
        ('C#CC(317)', 'C=CC(38)'),
    ],
    bathGas = {
        'Ar': 0.25,
        'He': 0.25,
        'Ne': 0.25,
        'N2': 0.25,
    },
)

pressureDependence(
    label = 'PDepNetwork #163',
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


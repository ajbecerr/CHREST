species(
    label = 'C[C]CC(167)',
    structure = SMILES('C[C]CC'),
    E0 = (284.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0351742,'amu*angstrom^2'), symmetry=1, barrier=(10.8251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.53834,'amu*angstrom^2'), symmetry=1, barrier=(58.3614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.53917,'amu*angstrom^2'), symmetry=1, barrier=(58.3804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4418.31,'J/mol'), sigma=(4.687e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23552,0.0362128,-1.62669e-05,2.7719e-09,-6.63706e-14,34266.3,25.257], Tmin=(100,'K'), Tmax=(1802.4,'K')), NASAPolynomial(coeffs=[10.5023,0.0226356,-8.93642e-06,1.5285e-09,-9.75175e-14,30511.7,-21.6443], Tmin=(1802.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=CCC(33)',
    structure = SMILES('C=CCC'),
    E0 = (-16.5218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,385.437],'cm^-1')),
        HinderedRotor(inertia=(0.105994,'amu*angstrom^2'), symmetry=1, barrier=(11.0704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10982,'amu*angstrom^2'), symmetry=1, barrier=(11.1095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64175,0.02037,3.05815e-05,-4.85565e-08,1.85234e-11,-1929.81,14.9537], Tmin=(100,'K'), Tmax=(991.791,'K')), NASAPolynomial(coeffs=[7.8083,0.0229242,-8.65891e-06,1.60051e-09,-1.13878e-13,-4105.08,-15.7292], Tmin=(991.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.5218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=CC(320)',
    structure = SMILES('CC=CC'),
    E0 = (-29.9019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.471166,'amu*angstrom^2'), symmetry=1, barrier=(10.833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471389,'amu*angstrom^2'), symmetry=1, barrier=(10.8382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61747,0.024342,1.17373e-05,-2.39959e-08,8.68014e-12,-3541.51,13.2442], Tmin=(100,'K'), Tmax=(1081.86,'K')), NASAPolynomial(coeffs=[5.99489,0.0261388,-1.05589e-05,1.94788e-09,-1.35284e-13,-5108.22,-7.179], Tmin=(1081.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.9019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = '[CH2]C[CH]C(34)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1553.59],'cm^-1')),
        HinderedRotor(inertia=(0.00260962,'amu*angstrom^2'), symmetry=1, barrier=(4.49698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195019,'amu*angstrom^2'), symmetry=1, barrier=(4.48387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19415,'amu*angstrom^2'), symmetry=1, barrier=(4.46389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2896.88,'J/mol'), sigma=(5.41451,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.49 K, Pc=41.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59222,0.0286124,-6.18623e-06,-3.09363e-09,1.21718e-12,30768.7,19.1945], Tmin=(100,'K'), Tmax=(1492.35,'K')), NASAPolynomial(coeffs=[6.22715,0.0257402,-1.02052e-05,1.78672e-09,-1.17177e-13,28918.7,-2.36048], Tmin=(1492.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2]C=C(53)',
    structure = SMILES('[CH2]C=C'),
    E0 = (157.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.570286,'amu*angstrom^2'), symmetry=1, barrier=(32.8574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31926,0.00566541,4.27429e-05,-5.78804e-08,2.21687e-11,18990.6,9.19662], Tmin=(100,'K'), Tmax=(952.01,'K')), NASAPolynomial(coeffs=[7.5573,0.0114809,-3.63938e-06,6.63549e-10,-4.95289e-14,17113.2,-16.6632], Tmin=(952.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(50)',
    structure = SMILES('[CH]=C'),
    E0 = (289.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,826.013,826.013,3240.27],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90668,-0.00406215,3.8677e-05,-4.62963e-08,1.72894e-11,34797.2,6.09797], Tmin=(100,'K'), Tmax=(931.97,'K')), NASAPolynomial(coeffs=[5.44803,0.00498345,-1.08814e-06,1.79821e-10,-1.45083e-14,33829.7,-4.87845], Tmin=(931.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C[CH2](7)',
    structure = SMILES('C[CH2]'),
    E0 = (108.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,474.393,1048.54,2320.41,2320.68,2321.07],'cm^-1')),
        HinderedRotor(inertia=(0.000750666,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82185,-0.00343378,5.09264e-05,-6.20221e-08,2.37078e-11,13066,7.61638], Tmin=(100,'K'), Tmax=(900.31,'K')), NASAPolynomial(coeffs=[5.15617,0.00943129,-1.8195e-06,2.21205e-10,-1.43489e-14,12064.1,-2.91078], Tmin=(900.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ)"""),
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
    label = '[CH]=CC(69)',
    structure = SMILES('[CH]=CC'),
    E0 = (253.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,180,1112.61,1112.63,3990.1],'cm^-1')),
        HinderedRotor(inertia=(0.181748,'amu*angstrom^2'), symmetry=1, barrier=(4.17875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.23408,0.0118208,1.70307e-05,-2.64368e-08,9.91228e-12,30487.3,10.3183], Tmin=(100,'K'), Tmax=(997.874,'K')), NASAPolynomial(coeffs=[5.6647,0.0144326,-5.46739e-06,1.00158e-09,-7.04862e-14,29387.1,-4.48505], Tmin=(997.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C=C[CH]C(164)',
    structure = SMILES('[CH2]C=CC'),
    E0 = (121.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.46744,'amu*angstrom^2'), symmetry=1, barrier=(10.7474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247286,'amu*angstrom^2'), symmetry=1, barrier=(24.7829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64215,0.0215952,2.09692e-05,-3.79221e-08,1.4785e-11,14680.9,14.1266], Tmin=(100,'K'), Tmax=(998.554,'K')), NASAPolynomial(coeffs=[7.82797,0.0208401,-7.96758e-06,1.47339e-09,-1.04526e-13,12647.3,-15.8818], Tmin=(998.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=CC(322)',
    structure = SMILES('C[C]=CC'),
    E0 = (207.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.423109,'amu*angstrom^2'), symmetry=1, barrier=(9.72812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423367,'amu*angstrom^2'), symmetry=1, barrier=(9.73403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62895,0.0282209,-9.66322e-06,-2.59435e-10,5.17284e-13,25060.3,14.3838], Tmin=(100,'K'), Tmax=(1614.54,'K')), NASAPolynomial(coeffs=[7.74501,0.0213061,-8.59049e-06,1.50735e-09,-9.84511e-14,22657.5,-15.0741], Tmin=(1614.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S)"""),
)

species(
    label = 'C[C]C(70)',
    structure = SMILES('C[C]C'),
    E0 = (308.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(2.2015,'amu*angstrom^2'), symmetry=1, barrier=(50.6167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20137,'amu*angstrom^2'), symmetry=1, barrier=(50.6137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1971.36,'J/mol'), sigma=(5.118e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92363,0.0210869,-1.82391e-06,-5.35956e-09,1.86097e-12,37101.5,19.8374], Tmin=(100,'K'), Tmax=(1330.24,'K')), NASAPolynomial(coeffs=[5.20765,0.0205522,-8.36243e-06,1.49635e-09,-1.00129e-13,35933.5,6.06042], Tmin=(1330.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH2]CC[CH2](23)',
    structure = SMILES('[CH2]CC[CH2]'),
    E0 = (266.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1262.28],'cm^-1')),
        HinderedRotor(inertia=(0.0854236,'amu*angstrom^2'), symmetry=1, barrier=(1.96406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854285,'amu*angstrom^2'), symmetry=1, barrier=(1.96417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854203,'amu*angstrom^2'), symmetry=1, barrier=(1.96398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2896.88,'J/mol'), sigma=(5.41451,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.49 K, Pc=41.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37803,0.0299378,-7.71059e-07,-1.37453e-08,5.72553e-12,32078.3,18.3087], Tmin=(100,'K'), Tmax=(1093.21,'K')), NASAPolynomial(coeffs=[7.15009,0.0244836,-9.76144e-06,1.78352e-09,-1.23052e-13,30317.5,-8.42047], Tmin=(1093.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=C(54)',
    structure = SMILES('[CH2]CC=C'),
    E0 = (188.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.067984,'amu*angstrom^2'), symmetry=1, barrier=(1.56309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223225,'amu*angstrom^2'), symmetry=1, barrier=(13.4708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6806,0.0210828,2.02115e-05,-3.64233e-08,1.4144e-11,22752.8,16.6009], Tmin=(100,'K'), Tmax=(1000.95,'K')), NASAPolynomial(coeffs=[7.5948,0.0206425,-7.89783e-06,1.45964e-09,-1.03413e-13,20807.3,-11.916], Tmin=(1000.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]CC(165)',
    structure = SMILES('C=[C]CC'),
    E0 = (221.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,350.561],'cm^-1')),
        HinderedRotor(inertia=(0.10577,'amu*angstrom^2'), symmetry=1, barrier=(9.20007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105602,'amu*angstrom^2'), symmetry=1, barrier=(9.20022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59444,0.0251166,5.39503e-06,-1.89578e-08,7.5085e-12,26674.2,15.5985], Tmin=(100,'K'), Tmax=(1056.58,'K')), NASAPolynomial(coeffs=[6.83079,0.0220763,-8.74093e-06,1.60417e-09,-1.11497e-13,25053.5,-8.50584], Tmin=(1056.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCC(166)',
    structure = SMILES('[CH]=CCC'),
    E0 = (230.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332475,'amu*angstrom^2'), symmetry=1, barrier=(11.7968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332229,'amu*angstrom^2'), symmetry=1, barrier=(11.7984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56532,0.0239246,1.38114e-05,-3.10267e-08,1.25443e-11,27790.1,15.6311], Tmin=(100,'K'), Tmax=(999.922,'K')), NASAPolynomial(coeffs=[8.06559,0.0200625,-7.60831e-06,1.39784e-09,-9.86868e-14,25783.2,-15.4389], Tmin=(999.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]CCC(168)',
    structure = SMILES('[CH]CCC'),
    E0 = (306.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,205.113,869.504,869.755,3356.71],'cm^-1')),
        HinderedRotor(inertia=(0.214808,'amu*angstrom^2'), symmetry=1, barrier=(6.42019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00400561,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215607,'amu*angstrom^2'), symmetry=1, barrier=(6.40955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24416,0.0351852,-1.67673e-05,3.10101e-09,-9.37033e-14,36965.1,16.0694], Tmin=(100,'K'), Tmax=(1696.97,'K')), NASAPolynomial(coeffs=[10.2265,0.021086,-8.47349e-06,1.48049e-09,-9.62463e-14,33576.9,-28.663], Tmin=(1696.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH2][CH2](18)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (313.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1608.02,1609.59,1610.81,1612.2,1612.86],'cm^-1')),
        HinderedRotor(inertia=(0.00575811,'amu*angstrom^2'), symmetry=1, barrier=(10.593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86546,-0.00370153,4.71859e-05,-6.15022e-08,2.54427e-11,37752.3,8.63154], Tmin=(100,'K'), Tmax=(845.133,'K')), NASAPolynomial(coeffs=[5.8996,0.00441368,1.29131e-06,-4.57987e-10,3.67865e-14,36774.8,-4.58851], Tmin=(845.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ)"""),
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
    label = 'CH2(T)(15)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.98,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[CH]C[CH2](312)',
    structure = SMILES('[CH]C[CH2]'),
    E0 = (532.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,1214.34,1214.82],'cm^-1')),
        HinderedRotor(inertia=(0.141587,'amu*angstrom^2'), symmetry=1, barrier=(3.25536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00310575,'amu*angstrom^2'), symmetry=1, barrier=(3.25017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04384,0.0162269,8.19929e-06,-1.95239e-08,7.97721e-12,64136.3,13.5821], Tmin=(100,'K'), Tmax=(994.051,'K')), NASAPolynomial(coeffs=[6.51415,0.0136529,-5.10446e-06,9.25401e-10,-6.47359e-14,62883.6,-5.97006], Tmin=(994.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C[C]C(498)',
    structure = SMILES('[CH2]C[C]C'),
    E0 = (509.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,180,870.619],'cm^-1')),
        HinderedRotor(inertia=(1.17608,'amu*angstrom^2'), symmetry=1, barrier=(27.0404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00100424,'amu*angstrom^2'), symmetry=1, barrier=(11.4022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80985,'amu*angstrom^2'), symmetry=1, barrier=(41.6119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28697,0.0322908,-1.00184e-05,-5.96281e-09,3.57489e-12,61303.7,17.4467], Tmin=(100,'K'), Tmax=(1088.51,'K')), NASAPolynomial(coeffs=[8.15768,0.0208543,-8.22742e-06,1.4956e-09,-1.02998e-13,59425.1,-14.1325], Tmin=(1088.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C[CH]C(499)',
    structure = SMILES('[CH]C[CH]C'),
    E0 = (498.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,210.638,1385.54,1385.54,1385.54],'cm^-1')),
        HinderedRotor(inertia=(0.17227,'amu*angstrom^2'), symmetry=1, barrier=(5.42427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00398198,'amu*angstrom^2'), symmetry=1, barrier=(5.42452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172279,'amu*angstrom^2'), symmetry=1, barrier=(5.42421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54331,0.0304752,-1.37918e-05,2.71973e-09,-1.80229e-13,59992.3,18.587], Tmin=(100,'K'), Tmax=(2031.72,'K')), NASAPolynomial(coeffs=[10.3933,0.0179392,-6.69153e-06,1.09704e-09,-6.75683e-14,56200.1,-26.3319], Tmin=(2031.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC1CC1(500)',
    structure = SMILES('CC1CC1'),
    E0 = (9.25087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02553,0.00382842,8.75459e-05,-1.14176e-07,4.3452e-11,1163.75,12.3752], Tmin=(100,'K'), Tmax=(948.048,'K')), NASAPolynomial(coeffs=[10.4141,0.0175357,-5.15238e-06,9.44366e-10,-7.28488e-14,-2254.14,-33.5135], Tmin=(948.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.25087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane)"""),
)

species(
    label = 'C=C(17)',
    structure = SMILES('C=C'),
    E0 = (42.1493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97974,-0.0075756,5.52973e-05,-6.36222e-08,2.31767e-11,5077.46,4.04623], Tmin=(100,'K'), Tmax=(940.446,'K')), NASAPolynomial(coeffs=[5.20299,0.00782443,-2.12684e-06,3.79691e-10,-2.94671e-14,3936.28,-6.6241], Tmin=(940.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.1493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2][CH][CH]C(501)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (449.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.000569455,'amu*angstrom^2'), symmetry=1, barrier=(2.14719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0174147,'amu*angstrom^2'), symmetry=1, barrier=(65.6618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000569937,'amu*angstrom^2'), symmetry=1, barrier=(2.14891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44689,0.0215254,-1.95614e-06,-2.95668e-09,7.24118e-13,54114,18.8308], Tmin=(100,'K'), Tmax=(2043.31,'K')), NASAPolynomial(coeffs=[8.78949,0.0205745,-8.23772e-06,1.37004e-09,-8.38825e-14,49945.9,-15.6186], Tmin=(2043.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH2](85)',
    structure = SMILES('[CH2][CH]C[CH2]'),
    E0 = (460.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2622.36],'cm^-1')),
        HinderedRotor(inertia=(0.0897296,'amu*angstrom^2'), symmetry=1, barrier=(2.06306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00870163,'amu*angstrom^2'), symmetry=1, barrier=(42.4425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0245193,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81948,0.0272005,-9.56974e-06,5.91729e-10,1.83286e-13,55442.9,20.1592], Tmin=(100,'K'), Tmax=(1946.59,'K')), NASAPolynomial(coeffs=[9.19037,0.0193116,-7.49963e-06,1.25572e-09,-7.83181e-14,51976.9,-17.352], Tmin=(1946.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC(35)',
    structure = SMILES('[CH2][CH]CC'),
    E0 = (255.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1553.59],'cm^-1')),
        HinderedRotor(inertia=(0.00260962,'amu*angstrom^2'), symmetry=1, barrier=(4.49698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195019,'amu*angstrom^2'), symmetry=1, barrier=(4.48387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19415,'amu*angstrom^2'), symmetry=1, barrier=(4.46389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59222,0.0286124,-6.18623e-06,-3.09363e-09,1.21718e-12,30768.7,19.1945], Tmin=(100,'K'), Tmax=(1492.35,'K')), NASAPolynomial(coeffs=[6.22715,0.0257402,-1.02052e-05,1.78672e-09,-1.17177e-13,28918.7,-2.36048], Tmin=(1492.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C(321)',
    structure = SMILES('C[CH][CH]C'),
    E0 = (244.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.00187131,'amu*angstrom^2'), symmetry=1, barrier=(6.82844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180225,'amu*angstrom^2'), symmetry=1, barrier=(65.775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00187326,'amu*angstrom^2'), symmetry=1, barrier=(6.82942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1236,0.0238473,-9.09881e-07,-4.52003e-09,1.14172e-12,29444.7,17.5352], Tmin=(100,'K'), Tmax=(1813.47,'K')), NASAPolynomial(coeffs=[6.29313,0.0260792,-1.03848e-05,1.76766e-09,-1.11708e-13,26778.5,-3.82402], Tmin=(1813.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC)"""),
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
    E0 = (425.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (278.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (388.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (333.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (419.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (321.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (727.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (371.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (425.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (277.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (329.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (293.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (397.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (333.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (400.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (433.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (442.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (343.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (722.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (717.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (725.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (720.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (710.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (262.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (344.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (403.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (387.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (393.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (661.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (672.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (412.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (396.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (418.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(2)', 'C=CC(38)'],
    products = ['CC=CC(320)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 4 CH2 + CH3CHCH2_r1 <=> CH3CHCHCH3 in 1,2_Insertion_carbene/training
This reaction matched rate rule [carbene;Cd_pri]
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C[CH]C(34)'],
    products = ['CC=CC(320)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH3](12)', '[CH]=CC(69)'],
    products = ['CC=CC(320)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.95927e+08,'m^3/(mol*s)'), n=-0.639165, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0851650107975, var=0.528125578779, Tref=1000.0, N=11, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-3R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-3R!H-R
    Total Standard Deviation in ln(k): 1.67086850823
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-3R!H-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(9)', 'C=C[CH]C(164)'],
    products = ['CC=CC(320)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(9)', 'C[C]=CC(322)'],
    products = ['CC=CC(320)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[C]CC(167)'],
    products = ['CC=CC(320)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(2)', 'C[C]C(70)'],
    products = ['C[C]CC(167)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.62043e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[C]CC(167)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.85495e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(S)(2)', 'C=CC(38)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.00284,'kcal/mol'), T0=(1,'K'), comment="""Matched reaction 6 CH2 + CH3CHCH2_r3 <=> CH3CH2CHCH2 in 1,2_Insertion_carbene/training
This reaction matched rate rule [carbene;C_pri/Cd]
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C[CH]C(34)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC[CH2](23)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH3](12)', '[CH2]C=C(53)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.04e+14,'cm^3/(mol*s)'), n=-0.32, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Matched reaction 127 C3H5 + CH3 <=> C4H8-2 in R_Recombination/training
This reaction matched rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-3R!H-R_N-Sp-3R!H=2R_Sp-4R!H=3R!H]
family: R_Recombination
Ea raised from -0.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(50)', 'C[CH2](7)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.08049e+08,'m^3/(mol*s)'), n=-0.685207, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.691788023649, var=1.05747121089, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Sp-3R!H-2R_3R!H->C_2R->C_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Sp-3R!H-2R_3R!H->C_2R->C_Ext-1C-R
    Total Standard Deviation in ln(k): 3.79969848986
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Sp-3R!H-2R_3R!H->C_2R->C_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(9)', 'C=C[CH]C(164)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Matched reaction 122 H + C4H7 <=> C4H8 in R_Recombination/training
This reaction matched rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-3R!H-R]
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(9)', '[CH2]CC=C(54)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.17499e+07,'m^3/(mol*s)'), n=0.115342, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(9)', 'C=[C]CC(165)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(9)', '[CH]=CCC(166)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]CCC(168)'],
    products = ['C=CCC(33)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH2](18)', '[CH]C(16)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.45208e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(T)(15)', '[CH2][CH]C(39)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH3](12)', '[CH]C[CH2](312)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(9)', '[CH2]C[C]C(498)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(9)', '[CH]C[CH]C(499)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C[CH]C(34)'],
    products = ['CC1CC1(500)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(9)', 'C=C[CH]C(164)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(170.395,'m^3/(mol*s)'), n=1.5621, Ea=(11.2886,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(9)', '[CH2]CC=C(54)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(T)(15)', 'C=CC(38)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06e+08,'cm^3/(mol*s)'), n=1.58, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-HH_Cds-Cs\H3/H;Y_1centerbirad] for rate rule [Cds-HH_Cds-Cs\H3/H;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -8.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C(17)', '[CH]C(16)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.84e+07,'cm^3/(mol*s)'), n=1.55, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(9)', '[CH2][CH][CH]C(501)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(9)', '[CH2][CH]C[CH2](85)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]CC(35)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C[CH]C(34)'],
    products = ['C[CH][CH]C(321)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC[CH2](23)'],
    products = ['[CH2]C[CH]C(34)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.436e+06,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #88',
    isomers = [
        'C[C]CC(167)',
        'C=CCC(33)',
        'CC=CC(320)',
        '[CH2]C[CH]C(34)',
    ],
    reactants = [
        ('[CH3](12)', '[CH2]C=C(53)'),
        ('[CH]=C(50)', 'C[CH2](7)'),
    ],
    bathGas = {
        'Ar': 0.25,
        'He': 0.25,
        'Ne': 0.25,
        'N2': 0.25,
    },
)

pressureDependence(
    label = 'PDepNetwork #88',
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


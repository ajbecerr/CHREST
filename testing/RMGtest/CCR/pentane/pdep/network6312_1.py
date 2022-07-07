species(
    label = '[CH2]C=CCC([CH2])([O])[CH][O](24788)',
    structure = SMILES('[CH2]C=CCC([CH2])([O])[CH][O]'),
    E0 = (548.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,180,180,180,180,774.431,1450.91,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154615,'amu*angstrom^2'), symmetry=1, barrier=(3.5549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154615,'amu*angstrom^2'), symmetry=1, barrier=(3.5549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154615,'amu*angstrom^2'), symmetry=1, barrier=(3.5549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154615,'amu*angstrom^2'), symmetry=1, barrier=(3.5549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154615,'amu*angstrom^2'), symmetry=1, barrier=(3.5549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.919123,0.116695,-0.00017077,1.40994e-07,-4.68613e-11,66083.6,37.4512], Tmin=(100,'K'), Tmax=(806.493,'K')), NASAPolynomial(coeffs=[11.8123,0.0449,-2.11478e-05,4.01303e-09,-2.7647e-13,64311.4,-19.4882], Tmin=(806.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4007.73,'J/mol'), sigma=(6.4029,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=626.00 K, Pc=34.64 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
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
    label = '[CH2]C=CCC([CH2])([O])C=O(24800)',
    structure = SMILES('[CH2]C=CCC([CH2])([O])C=O'),
    E0 = (242.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1631.59,2875.34,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150681,'amu*angstrom^2'), symmetry=1, barrier=(3.46445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150681,'amu*angstrom^2'), symmetry=1, barrier=(3.46445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150681,'amu*angstrom^2'), symmetry=1, barrier=(3.46445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150681,'amu*angstrom^2'), symmetry=1, barrier=(3.46445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150681,'amu*angstrom^2'), symmetry=1, barrier=(3.46445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.474525,0.100597,-0.00011463,7.04476e-08,-1.74539e-11,29286.5,36.5869], Tmin=(100,'K'), Tmax=(979.201,'K')), NASAPolynomial(coeffs=[15.3901,0.0357909,-1.53557e-05,2.859e-09,-1.97833e-13,26179.5,-39.6176], Tmin=(979.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CC(C)(C=O)OJ) + radical(CJC(C)(C=O)O) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]CC[C]([O])[CH][O](25680)',
    structure = SMILES('[CH2]C=CCC[C]([O])[CH][O]'),
    E0 = (532.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4567.69,'J/mol'), sigma=(7.66872,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=713.46 K, Pc=22.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.598696,0.11052,-0.000162226,1.38799e-07,-4.78666e-11,64151.6,37.8641], Tmin=(100,'K'), Tmax=(812.897,'K')), NASAPolynomial(coeffs=[9.29829,0.0484453,-2.30021e-05,4.37904e-09,-3.02092e-13,62984.5,-5.11487], Tmin=(812.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]([O])C([O])C[CH]C=C(25655)',
    structure = SMILES('[CH2]C=CCC([O])[C]([CH2])[O]'),
    E0 = (555.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.476475,0.10094,-0.000111622,6.62288e-08,-1.59137e-11,66957.1,37.1228], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[15.5506,0.037166,-1.6459e-05,3.11734e-09,-2.18006e-13,63734.9,-40.2822], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CC[C]([CH2])[CH][O](27537)',
    structure = SMILES('[CH2]C=CC[C]([CH2])[CH][O]'),
    E0 = (646.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656739,0.0777146,-7.57729e-05,4.53514e-08,-1.17965e-11,77834.5,33.4795], Tmin=(100,'K'), Tmax=(904.099,'K')), NASAPolynomial(coeffs=[8.37046,0.0435855,-1.91467e-05,3.59469e-09,-2.49508e-13,76439.7,-2.95676], Tmin=(904.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CCsJOH) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]CC([O])=C[O](11358)',
    structure = SMILES('[CH2][CH][CH]CC([O])=C[O]'),
    E0 = (385.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,201.777,802.666,1070.22,1342.22,1614.22],'cm^-1')),
        HinderedRotor(inertia=(0.152936,'amu*angstrom^2'), symmetry=1, barrier=(3.57907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152936,'amu*angstrom^2'), symmetry=1, barrier=(3.57907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152936,'amu*angstrom^2'), symmetry=1, barrier=(3.57907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152936,'amu*angstrom^2'), symmetry=1, barrier=(3.57907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.668674,0.0738126,-7.83344e-05,4.67081e-08,-1.13642e-11,46542.8,36.4427], Tmin=(100,'K'), Tmax=(993.069,'K')), NASAPolynomial(coeffs=[11.5811,0.0298592,-1.19454e-05,2.1406e-09,-1.44737e-13,44375.4,-16.1276], Tmin=(993.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH][O](1548)',
    structure = SMILES('[CH][O]'),
    E0 = (431.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([806.798,806.798,2696.84],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86882,-0.000695042,1.53508e-05,-2.00576e-08,7.65929e-12,51865.1,6.76618], Tmin=(100,'K'), Tmax=(973.285,'K')), NASAPolynomial(coeffs=[6.03251,-0.000301449,4.32946e-07,-3.67045e-11,-1.25373e-15,51004.1,-5.87326], Tmin=(973.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2439,0.0656708,-7.98273e-05,6.56523e-08,-2.26796e-11,54615.6,33.1053], Tmin=(100,'K'), Tmax=(830.908,'K')), NASAPolynomial(coeffs=[4.97977,0.0392941,-1.70605e-05,3.13676e-09,-2.12872e-13,54284.4,17.5173], Tmin=(830.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CCC([CH2])([O])[CH][O](27538)',
    structure = SMILES('[CH]=CCC([CH2])([O])[CH][O]'),
    E0 = (679.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,364.295,796.24,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.41755,0.108013,-0.000178594,1.57479e-07,-5.37744e-11,81893.8,33.9648], Tmin=(100,'K'), Tmax=(853.569,'K')), NASAPolynomial(coeffs=[10.1224,0.0376712,-1.8167e-05,3.42645e-09,-2.32832e-13,80857.6,-10.7451], Tmin=(853.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Cds_P)"""),
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
    label = '[CH]C([O])([CH][O])CC=C[CH2](27539)',
    structure = SMILES('[CH]C([O])([CH][O])CC=C[CH2]'),
    E0 = (784.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.38886,0.105315,-0.000136344,9.26935e-08,-2.27028e-11,94472.7,36.7513], Tmin=(100,'K'), Tmax=(651.9,'K')), NASAPolynomial(coeffs=[12.0912,0.041665,-1.96337e-05,3.75764e-09,-2.61657e-13,92570.9,-20.2248], Tmin=(651.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(784.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC([CH2])([O])[C][O](27540)',
    structure = SMILES('[CH2]C=CCC([CH2])([O])[C][O]'),
    E0 = (828.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.782338,0.112825,-0.000160278,1.25126e-07,-3.94614e-11,99838.9,36.1514], Tmin=(100,'K'), Tmax=(774.247,'K')), NASAPolynomial(coeffs=[13.3054,0.0400434,-1.92757e-05,3.7168e-09,-2.59481e-13,97657.4,-28.2096], Tmin=(774.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(828.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C[CH]CC([CH2])([O])[CH][O](25909)',
    structure = SMILES('[CH]C=CCC([CH2])([O])[CH][O]'),
    E0 = (767.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.943502,0.119432,-0.00018018,1.55571e-07,-5.345e-11,92444.2,38.365], Tmin=(100,'K'), Tmax=(834.251,'K')), NASAPolynomial(coeffs=[9.52764,0.0509236,-2.4093e-05,4.54119e-09,-3.10416e-13,91334,-6.4374], Tmin=(834.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC1([CH][O])CO1(27541)',
    structure = SMILES('[CH2]C=CCC1([CH][O])CO1'),
    E0 = (296.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.583215,0.0980223,-0.000109303,6.82853e-08,-1.70044e-11,35849.3,33.3804], Tmin=(100,'K'), Tmax=(1018.32,'K')), NASAPolynomial(coeffs=[15.1758,0.0337522,-1.11435e-05,1.73934e-09,-1.06533e-13,32762.6,-42.3307], Tmin=(1018.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1(C[CH]C=C)OC1[O](24062)',
    structure = SMILES('[CH2]C=CCC1([CH2])OC1[O]'),
    E0 = (286.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4308.76,'J/mol'), sigma=(7.29087,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=673.02 K, Pc=25.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.426131,0.094847,-9.41817e-05,4.63066e-08,-7.10337e-12,34618.1,32.7931], Tmin=(100,'K'), Tmax=(850.451,'K')), NASAPolynomial(coeffs=[16.0196,0.0325317,-1.07912e-05,1.72522e-09,-1.0912e-13,31277.1,-47.0809], Tmin=(850.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC1([O])CC1[O](27542)',
    structure = SMILES('[CH2]C=CCC1([O])CC1[O]'),
    E0 = (310.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.175748,0.0758608,-3.01814e-05,-2.10245e-08,1.48927e-11,37489.8,32.8551], Tmin=(100,'K'), Tmax=(988.986,'K')), NASAPolynomial(coeffs=[20.8505,0.0266122,-9.7734e-06,1.81342e-09,-1.30976e-13,31580.4,-77.2013], Tmin=(988.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([CH][O])CC=CCO1(25668)',
    structure = SMILES('[CH2]C1([CH][O])CC=CCO1'),
    E0 = (258.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.097755,0.0892231,-8.78404e-05,4.91923e-08,-1.13242e-11,31209,28.2545], Tmin=(100,'K'), Tmax=(1044.31,'K')), NASAPolynomial(coeffs=[13.1609,0.0384376,-1.48927e-05,2.6229e-09,-1.75605e-13,28439.8,-36.2855], Tmin=(1044.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1([O])CC=CCC1(25682)',
    structure = SMILES('[O][CH]C1([O])CC=CCC1'),
    E0 = (239.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0105116,0.0759595,-5.14569e-05,1.72313e-08,-2.31591e-12,28969.3,30.4153], Tmin=(100,'K'), Tmax=(1731.57,'K')), NASAPolynomial(coeffs=[19.5249,0.0308801,-1.24058e-05,2.19618e-09,-1.45165e-13,22211.2,-74.4441], Tmin=(1731.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])CC=CCC1[O](25657)',
    structure = SMILES('[CH2]C1([O])CC=CCC1[O]'),
    E0 = (263.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0025694,0.0721775,-2.25719e-05,-2.47253e-08,1.49714e-11,31828,28.6269], Tmin=(100,'K'), Tmax=(1008.58,'K')), NASAPolynomial(coeffs=[19.3387,0.0301704,-1.1704e-05,2.20268e-09,-1.58617e-13,26161.7,-73.5983], Tmin=(1008.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(O)([CH][O])CC=C=C(25672)',
    structure = SMILES('[CH2]C(O)([CH][O])CC=C=C'),
    E0 = (344.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14645,0.123403,-0.000192675,1.6357e-07,-5.45301e-11,41577.3,36.9772], Tmin=(100,'K'), Tmax=(845.667,'K')), NASAPolynomial(coeffs=[12.2224,0.0431532,-2.01517e-05,3.77217e-09,-2.56222e-13,39924.6,-21.6813], Tmin=(845.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH)"""),
)

species(
    label = 'C=C=CCC(C)([O])[CH][O](25683)',
    structure = SMILES('C=C=CCC(C)([O])[CH][O]'),
    E0 = (361.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.559633,0.108138,-0.000151254,1.23881e-07,-4.15234e-11,43591.1,35.0961], Tmin=(100,'K'), Tmax=(787.056,'K')), NASAPolynomial(coeffs=[10.341,0.0464844,-2.18323e-05,4.15926e-09,-2.88008e-13,42068.9,-13.6522], Tmin=(787.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])(C[O])CC=C=C(25662)',
    structure = SMILES('[CH2]C([O])(C[O])CC=C=C'),
    E0 = (392.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.630109,0.110344,-0.000158159,1.3205e-07,-4.48364e-11,47407.8,35.6371], Tmin=(100,'K'), Tmax=(796.365,'K')), NASAPolynomial(coeffs=[10.1622,0.0471929,-2.23636e-05,4.26837e-09,-2.95479e-13,45972.5,-12.1916], Tmin=(796.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])([O])[CH][O](27543)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])([O])[CH][O]'),
    E0 = (883.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 8,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.898841,0.123193,-0.000208302,1.91695e-07,-6.75317e-11,106476,44.0771], Tmin=(100,'K'), Tmax=(864.291,'K')), NASAPolynomial(coeffs=[7.37766,0.0510648,-2.44178e-05,4.57712e-09,-3.09367e-13,106309,12.6634], Tmin=(864.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(883.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCJCO) + radical(RCCJC) + radical(CJC(C)2O) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([O])([CH][O])CC1[CH]C1(27544)',
    structure = SMILES('[CH2]C([O])([CH][O])CC1[CH]C1'),
    E0 = (663.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.285935,0.10128,-0.000129026,9.71903e-08,-3.051e-11,79987.4,36.4947], Tmin=(100,'K'), Tmax=(769.8,'K')), NASAPolynomial(coeffs=[10.3242,0.0461456,-2.15869e-05,4.13987e-09,-2.89255e-13,78354,-11.9164], Tmin=(769.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(663.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(cyclopropane) + radical(CJC(C)2O) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])([CH][O])O1(27401)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])([CH][O])O1'),
    E0 = (553.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.746958,0.102048,-0.000120612,7.83178e-08,-1.99775e-11,66794.8,32.7817], Tmin=(100,'K'), Tmax=(1036.24,'K')), NASAPolynomial(coeffs=[15.9644,0.0315329,-9.84403e-06,1.46087e-09,-8.56281e-14,63653.9,-46.8796], Tmin=(1036.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCOJ) + radical(CCJCO) + radical(CJC(C)OC) + radical(CCsJOH) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1[CH]CC([O])([CH][O])C1(27235)',
    structure = SMILES('[CH2]C1[CH]CC([O])([CH][O])C1'),
    E0 = (535.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.367466,0.0744258,-5.36789e-05,2.0726e-08,-3.3301e-12,64599.3,35.8033], Tmin=(100,'K'), Tmax=(1435.5,'K')), NASAPolynomial(coeffs=[13.7647,0.0370944,-1.46697e-05,2.60942e-09,-1.74968e-13,60753,-33.6736], Tmin=(1435.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(Cs_S) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])([O])C1[O](27509)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])([O])C1[O]'),
    E0 = (559.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0696668,0.073569,-3.31715e-05,-1.27122e-08,1.12139e-11,67471.4,35.0675], Tmin=(100,'K'), Tmax=(980.471,'K')), NASAPolynomial(coeffs=[17.7527,0.0301576,-1.07102e-05,1.90094e-09,-1.32384e-13,62622.9,-56.9365], Tmin=(980.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(Cs_S) + radical(CJC(C)2O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])([O])[CH]O1(27331)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])([O])[CH]O1'),
    E0 = (536.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45279,0.101969,-0.000103818,5.4459e-08,-1.08643e-11,64795.5,33.5486], Tmin=(100,'K'), Tmax=(1397.92,'K')), NASAPolynomial(coeffs=[22.041,0.0218956,-4.11093e-06,3.34016e-10,-8.94824e-15,59482.5,-83.175], Tmin=(1397.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxane) + radical(CC(C)2OJ) + radical(CCJCO) + radical(CCsJOCs) + radical(CJC(C)2O) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])([CH][O])O1(25794)',
    structure = SMILES('[CH2][CH]C1CC([CH2])([CH][O])O1'),
    E0 = (627.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.292963,0.096378,-0.000100735,4.8835e-08,-4.73663e-12,75578.8,35.4206], Tmin=(100,'K'), Tmax=(742.089,'K')), NASAPolynomial(coeffs=[14.7315,0.0337482,-1.12417e-05,1.76817e-09,-1.0912e-13,72843.5,-35.9876], Tmin=(742.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(CCJCO) + radical(CJC(C)OC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([O])([CH][O])C1(27545)',
    structure = SMILES('[CH2][CH]C1CC([O])([CH][O])C1'),
    E0 = (617.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626615,0.0778308,-6.04757e-05,2.53554e-08,-4.56876e-12,74436.3,34.6827], Tmin=(100,'K'), Tmax=(1250.01,'K')), NASAPolynomial(coeffs=[11.1183,0.0442576,-2.0188e-05,3.86866e-09,-2.71432e-13,71813.4,-18.2746], Tmin=(1250.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(Cs_S) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])([O])C1[O](27546)',
    structure = SMILES('[CH2][CH]C1CC([CH2])([O])C1[O]'),
    E0 = (641.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.4816,0.0861817,-7.04666e-05,2.8953e-08,-4.73024e-12,77344,36.8765], Tmin=(100,'K'), Tmax=(1465.58,'K')), NASAPolynomial(coeffs=[20.7099,0.0283441,-1.12709e-05,2.02597e-09,-1.37026e-13,71132.4,-73.4607], Tmin=(1465.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(Cs_S) + radical(CJC(C)2O) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])([O])[CH]O1(27547)',
    structure = SMILES('[CH2][CH]C1CC([CH2])([O])[CH]O1'),
    E0 = (553.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39123,0.0953313,-9.02797e-05,4.40167e-08,-8.13068e-12,66789.1,36.3218], Tmin=(100,'K'), Tmax=(1536.64,'K')), NASAPolynomial(coeffs=[21.4007,0.0213579,-3.77519e-06,2.85205e-10,-6.89575e-15,61513.4,-77.8021], Tmin=(1536.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)2OJ) + radical(CCJCO) + radical(CCsJOCs) + radical(CJC(C)2O) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([O])([CH][O])CC=C=C(27548)',
    structure = SMILES('[CH2]C([O])([CH][O])CC=C=C'),
    E0 = (573.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180,180,919.997,1311.57,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160671,'amu*angstrom^2'), symmetry=1, barrier=(3.69414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160671,'amu*angstrom^2'), symmetry=1, barrier=(3.69414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160671,'amu*angstrom^2'), symmetry=1, barrier=(3.69414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160671,'amu*angstrom^2'), symmetry=1, barrier=(3.69414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.926948,0.119429,-0.000190619,1.65911e-07,-5.66002e-11,69100.8,36.2498], Tmin=(100,'K'), Tmax=(840.062,'K')), NASAPolynomial(coeffs=[10.9105,0.0437623,-2.10449e-05,3.988e-09,-2.72712e-13,67793,-14.7424], Tmin=(840.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O)"""),
)

species(
    label = 'C=C[CH]CC([O])=C[O](11260)',
    structure = SMILES('[CH2]C=CCC([O])=C[O]'),
    E0 = (60.1396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,371.32,371.337,371.345,371.374],'cm^-1')),
        HinderedRotor(inertia=(0.194783,'amu*angstrom^2'), symmetry=1, barrier=(19.0688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194838,'amu*angstrom^2'), symmetry=1, barrier=(19.0685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194864,'amu*angstrom^2'), symmetry=1, barrier=(19.0684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4496.12,'J/mol'), sigma=(7.19561,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=702.28 K, Pc=27.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315737,0.0648577,-1.99561e-05,-3.38401e-08,2.14782e-11,7380.62,30.9617], Tmin=(100,'K'), Tmax=(937.203,'K')), NASAPolynomial(coeffs=[21.8558,0.0136192,-3.08064e-06,4.86701e-10,-3.7325e-14,1555.93,-81.0942], Tmin=(937.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.1396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CCC(=C)[CH][O](22101)',
    structure = SMILES('[CH2]C=CCC([CH2])=C[O]'),
    E0 = (249.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,396.738,396.836,397.068],'cm^-1')),
        HinderedRotor(inertia=(0.217071,'amu*angstrom^2'), symmetry=1, barrier=(24.2167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216727,'amu*angstrom^2'), symmetry=1, barrier=(24.2192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216232,'amu*angstrom^2'), symmetry=1, barrier=(24.2136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217199,'amu*angstrom^2'), symmetry=1, barrier=(24.2163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117757,0.067472,-1.13368e-05,-4.28202e-08,2.37387e-11,30144.4,31.1511], Tmin=(100,'K'), Tmax=(957.795,'K')), NASAPolynomial(coeffs=[21.5968,0.0211315,-6.67161e-06,1.19989e-09,-8.87224e-14,24040.9,-81.9296], Tmin=(957.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]([O])[CH][O](26240)',
    structure = SMILES('[CH2][C]([O])[CH][O]'),
    E0 = (571.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,180,2035.01],'cm^-1')),
        HinderedRotor(inertia=(0.0866393,'amu*angstrom^2'), symmetry=1, barrier=(1.99201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086712,'amu*angstrom^2'), symmetry=1, barrier=(1.99368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60835,0.065059,-0.00013406,1.3209e-07,-4.73466e-11,68841.1,21.2809], Tmin=(100,'K'), Tmax=(890.641,'K')), NASAPolynomial(coeffs=[4.29645,0.0226555,-1.15626e-05,2.16119e-09,-1.43083e-13,69565.2,15.3771], Tmin=(890.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=C[CH2](16804)',
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
    label = '[CH2]C([CH2])([O])[CH][O](26275)',
    structure = SMILES('[CH2]C([CH2])([O])[CH][O]'),
    E0 = (564.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,180,180,180,280.802,1488.31,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.13832,'amu*angstrom^2'), symmetry=1, barrier=(3.18024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13832,'amu*angstrom^2'), symmetry=1, barrier=(3.18024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13832,'amu*angstrom^2'), symmetry=1, barrier=(3.18024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440529,0.0919539,-0.00017559,1.64977e-07,-5.75825e-11,67992.7,24.6642], Tmin=(100,'K'), Tmax=(884.555,'K')), NASAPolynomial(coeffs=[7.32782,0.0288663,-1.44403e-05,2.69762e-09,-1.79449e-13,68023.9,-0.654324], Tmin=(884.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])([O])[CH][O](27549)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])([O])[CH][O]'),
    E0 = (664.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.757828,0.115224,-0.000153912,1.01686e-07,-2.16877e-11,80137.9,34.1258], Tmin=(100,'K'), Tmax=(630.166,'K')), NASAPolynomial(coeffs=[13.5951,0.0415347,-1.99622e-05,3.8341e-09,-2.66791e-13,77983.1,-31.2348], Tmin=(630.166,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])([O])[CH][O](27550)',
    structure = SMILES('[CH2]C=[C]CC([CH2])([O])[CH][O]'),
    E0 = (785.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,776.092,1454.22,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00818,0.121994,-0.000198162,1.73817e-07,-5.93875e-11,94689.3,38.2418], Tmin=(100,'K'), Tmax=(847.606,'K')), NASAPolynomial(coeffs=[10.865,0.0439723,-2.1172e-05,4.00095e-09,-2.72656e-13,93466.5,-12.4164], Tmin=(847.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])([O])[CH][O](27551)',
    structure = SMILES('[CH2][C]=CCC([CH2])([O])[CH][O]'),
    E0 = (785.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,1685,370,180,180,180,180,776.092,1454.22,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154801,'amu*angstrom^2'), symmetry=1, barrier=(3.55918,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00818,0.121994,-0.000198162,1.73817e-07,-5.93875e-11,94689.3,38.2418], Tmin=(100,'K'), Tmax=(847.606,'K')), NASAPolynomial(coeffs=[10.865,0.0439723,-2.1172e-05,4.00095e-09,-2.72656e-13,93466.5,-12.4164], Tmin=(847.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])([CH][O])CC=[C]C(27552)',
    structure = SMILES('[CH2]C([O])([CH][O])CC=[C]C'),
    E0 = (634.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01271,0.124567,-0.000207131,1.87988e-07,-6.58949e-11,76465.9,37.9765], Tmin=(100,'K'), Tmax=(853.706,'K')), NASAPolynomial(coeffs=[8.60234,0.049973,-2.41567e-05,4.5664e-09,-3.10836e-13,75900.8,-0.584441], Tmin=(853.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(634.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])(O)[CH][O](27553)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])(O)[CH][O]'),
    E0 = (436.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29769,0.123848,-0.00017746,1.37799e-07,-4.28437e-11,52628.1,35.9535], Tmin=(100,'K'), Tmax=(787.94,'K')), NASAPolynomial(coeffs=[15.1138,0.0405242,-1.88156e-05,3.55442e-09,-2.44751e-13,50042.2,-39.3088], Tmin=(787.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJCO) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)([O])[CH][O](27554)',
    structure = SMILES('[CH2]C=C[CH]C(C)([O])[CH][O]'),
    E0 = (452.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.548187,0.106494,-0.000127816,8.5855e-08,-2.3743e-11,54634.9,33.4993], Tmin=(100,'K'), Tmax=(872.327,'K')), NASAPolynomial(coeffs=[13.0352,0.044205,-2.07036e-05,3.99153e-09,-2.80746e-13,52265.2,-30.1768], Tmin=(872.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])([O])C[O](27555)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])([O])C[O]'),
    E0 = (484.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.621185,0.108706,-0.000134627,9.37435e-08,-2.6871e-11,58451.7,34.0509], Tmin=(100,'K'), Tmax=(842.979,'K')), NASAPolynomial(coeffs=[12.7993,0.0450216,-2.13021e-05,4.1174e-09,-2.89665e-13,56189.2,-28.4024], Tmin=(842.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CJC(C)2O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([O])([CH][O])C[C]=CC(27556)',
    structure = SMILES('[CH2]C([O])([CH][O])C[C]=CC'),
    E0 = (634.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,180,180,180,180,1006.92,1231.03,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162862,'amu*angstrom^2'), symmetry=1, barrier=(3.74452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162862,'amu*angstrom^2'), symmetry=1, barrier=(3.74452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162862,'amu*angstrom^2'), symmetry=1, barrier=(3.74452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162862,'amu*angstrom^2'), symmetry=1, barrier=(3.74452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162862,'amu*angstrom^2'), symmetry=1, barrier=(3.74452,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01271,0.124567,-0.000207131,1.87988e-07,-6.58949e-11,76465.9,37.9765], Tmin=(100,'K'), Tmax=(853.706,'K')), NASAPolynomial(coeffs=[8.60234,0.049973,-2.41567e-05,4.5664e-09,-3.10836e-13,75900.8,-0.584441], Tmin=(853.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(634.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])([O])[CH]O(27557)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])([O])[CH]O'),
    E0 = (439.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.30558,0.120928,-0.000156638,1.05691e-07,-2.82819e-11,53017.8,35.1884], Tmin=(100,'K'), Tmax=(914.795,'K')), NASAPolynomial(coeffs=[18.3576,0.0349498,-1.56594e-05,2.95099e-09,-2.04772e-13,49420.2,-57.9242], Tmin=(914.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])(O)[CH][O](27558)',
    structure = SMILES('[CH2]C=[C]CC([CH2])(O)[CH][O]'),
    E0 = (556.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,180,180,239.773,916.391,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22833,0.125977,-0.000200259,1.71545e-07,-5.73569e-11,67165.9,38.9715], Tmin=(100,'K'), Tmax=(853.451,'K')), NASAPolynomial(coeffs=[12.1752,0.043366,-2.02804e-05,3.78551e-09,-2.56198e-13,65598.8,-19.3461], Tmin=(853.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(C)([O])[CH][O](27559)',
    structure = SMILES('[CH2]C=[C]CC(C)([O])[CH][O]'),
    E0 = (573.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,180,180,180,840.015,1389.46,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.658997,0.110937,-0.000159704,1.33091e-07,-4.49198e-11,69180.3,37.1519], Tmin=(100,'K'), Tmax=(805.587,'K')), NASAPolynomial(coeffs=[10.3545,0.0465881,-2.18954e-05,4.15667e-09,-2.86633e-13,67719.4,-11.6548], Tmin=(805.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])([O])C[O](27560)',
    structure = SMILES('[CH2]C=[C]CC([CH2])([O])C[O]'),
    E0 = (605.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,809.318,1420.8,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723846,0.113068,-0.00016631,1.40808e-07,-4.80088e-11,72996.8,37.6732], Tmin=(100,'K'), Tmax=(811.606,'K')), NASAPolynomial(coeffs=[10.1621,0.0473213,-2.24418e-05,4.26946e-09,-2.94418e-13,71628.2,-10.1188], Tmin=(811.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])([CH][O])[CH]C=CC(27561)',
    structure = SMILES('[CH2]C([O])([CH][O])[CH]C=CC'),
    E0 = (513.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14837,0.123308,-0.000187831,1.59526e-07,-5.40183e-11,61931,35.1909], Tmin=(100,'K'), Tmax=(817.258,'K')), NASAPolynomial(coeffs=[11.6877,0.046865,-2.25312e-05,4.29611e-09,-2.96057e-13,60287.7,-21.363], Tmin=(817.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CJC(C)2O) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])([O])[CH]O(27562)',
    structure = SMILES('[CH2]C=[C]CC([CH2])([O])[CH]O'),
    E0 = (560.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,180,180,180,265.418,881.327,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31758,0.124066,-0.000183245,1.44906e-07,-4.54382e-11,67559,38.4952], Tmin=(100,'K'), Tmax=(823.264,'K')), NASAPolynomial(coeffs=[15.4429,0.037757,-1.71072e-05,3.17859e-09,-2.15968e-13,64964.5,-38.1018], Tmin=(823.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])(O)[CH][O](27563)',
    structure = SMILES('[CH2][C]=CCC([CH2])(O)[CH][O]'),
    E0 = (556.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,1685,370,180,180,180,239.773,916.391,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151505,'amu*angstrom^2'), symmetry=1, barrier=(3.48339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22833,0.125977,-0.000200259,1.71545e-07,-5.73569e-11,67165.9,38.9715], Tmin=(100,'K'), Tmax=(853.451,'K')), NASAPolynomial(coeffs=[12.1752,0.043366,-2.02804e-05,3.78551e-09,-2.56198e-13,65598.8,-19.3461], Tmin=(853.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC(C)([O])[CH][O](27564)',
    structure = SMILES('[CH2][C]=CCC(C)([O])[CH][O]'),
    E0 = (573.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,180,180,180,180,840.015,1389.46,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156949,'amu*angstrom^2'), symmetry=1, barrier=(3.60857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.658997,0.110937,-0.000159704,1.33091e-07,-4.49198e-11,69180.3,37.1519], Tmin=(100,'K'), Tmax=(805.587,'K')), NASAPolynomial(coeffs=[10.3545,0.0465881,-2.18954e-05,4.15667e-09,-2.86633e-13,67719.4,-11.6548], Tmin=(805.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])([O])C[O](27565)',
    structure = SMILES('[CH2][C]=CCC([CH2])([O])C[O]'),
    E0 = (605.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180,180,180,809.318,1420.8,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155934,'amu*angstrom^2'), symmetry=1, barrier=(3.58524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723846,0.113068,-0.00016631,1.40808e-07,-4.80088e-11,72996.8,37.6732], Tmin=(100,'K'), Tmax=(811.606,'K')), NASAPolynomial(coeffs=[10.1621,0.0473213,-2.24418e-05,4.26946e-09,-2.94418e-13,71628.2,-10.1188], Tmin=(811.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])([O])[CH]O(27566)',
    structure = SMILES('[CH2][C]=CCC([CH2])([O])[CH]O'),
    E0 = (560.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180,180,265.418,881.327,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151883,'amu*angstrom^2'), symmetry=1, barrier=(3.49208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31758,0.124066,-0.000183245,1.44906e-07,-4.54382e-11,67559,38.4952], Tmin=(100,'K'), Tmax=(823.264,'K')), NASAPolynomial(coeffs=[15.4429,0.037757,-1.71072e-05,3.17859e-09,-2.15968e-13,64964.5,-38.1018], Tmin=(823.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)C([CH2])([O])[CH][O](24710)',
    structure = SMILES('[CH2]C(C=C)C([CH2])([O])[CH][O]'),
    E0 = (609.249,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,180,180,180,724.476,1508.31,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.1531,'amu*angstrom^2'), symmetry=1, barrier=(3.52006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1531,'amu*angstrom^2'), symmetry=1, barrier=(3.52006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1531,'amu*angstrom^2'), symmetry=1, barrier=(3.52006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1531,'amu*angstrom^2'), symmetry=1, barrier=(3.52006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1531,'amu*angstrom^2'), symmetry=1, barrier=(3.52006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4528.99,'J/mol'), sigma=(7.66898,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=707.42 K, Pc=22.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.990642,0.118493,-0.000177707,1.47023e-07,-4.79659e-11,73447.1,39.4016], Tmin=(100,'K'), Tmax=(856.748,'K')), NASAPolynomial(coeffs=[12.1251,0.0425497,-1.89935e-05,3.48421e-09,-2.3406e-13,71739.5,-18.6962], Tmin=(856.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Isobutyl)"""),
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
    label = '[CH]CC([CH2])([O])[CH][O](26278)',
    structure = SMILES('[CH]CC([CH2])([O])[CH][O]'),
    E0 = (780.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0130474,0.100007,-0.000175169,1.58379e-07,-5.44883e-11,93961.5,30.3951], Tmin=(100,'K'), Tmax=(869.273,'K')), NASAPolynomial(coeffs=[9.039,0.0331101,-1.61741e-05,3.03574e-09,-2.04436e-13,93341.5,-6.52178], Tmin=(869.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(780.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1([CH][O])CC(C=C)O1(25666)',
    structure = SMILES('[CH2]C1([CH][O])CC(C=C)O1'),
    E0 = (349.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.344328,0.0923065,-8.75972e-05,3.91858e-08,-4.37285e-12,42185.4,32.6834], Tmin=(100,'K'), Tmax=(851.709,'K')), NASAPolynomial(coeffs=[16.138,0.0316103,-1.01335e-05,1.58979e-09,-9.96549e-14,38771.6,-47.747], Tmin=(851.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = 'C=CC1CC([O])([CH][O])C1(27229)',
    structure = SMILES('C=CC1CC([O])([CH][O])C1'),
    E0 = (346.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.121585,0.0796886,-5.9108e-05,2.2137e-08,-3.33613e-12,41889.5,33.5549], Tmin=(100,'K'), Tmax=(1561.61,'K')), NASAPolynomial(coeffs=[18.8275,0.0311514,-1.24859e-05,2.23368e-09,-1.49802e-13,35971.2,-66.3097], Tmin=(1561.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])CC(C=C)C1[O](27500)',
    structure = SMILES('[CH2]C1([O])CC(C=C)C1[O]'),
    E0 = (370.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196355,0.076447,-3.13372e-05,-1.93132e-08,1.4073e-11,44751.4,32.0016], Tmin=(100,'K'), Tmax=(995.189,'K')), NASAPolynomial(coeffs=[20.7022,0.0274549,-1.02572e-05,1.91117e-09,-1.37847e-13,38858.2,-77.4307], Tmin=(995.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(O)([CH][O])C=CC=C(27567)',
    structure = SMILES('[CH2]C(O)([CH][O])C=CC=C'),
    E0 = (273.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.657533,0.105714,-0.000130216,8.66705e-08,-2.31192e-11,33028.3,37.2256], Tmin=(100,'K'), Tmax=(914.642,'K')), NASAPolynomial(coeffs=[15.3868,0.035545,-1.5137e-05,2.78874e-09,-1.91079e-13,30093.4,-38.7471], Tmin=(914.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=CC=CC(C)([O])[CH][O](27568)',
    structure = SMILES('C=CC=CC(C)([O])[CH][O]'),
    E0 = (288.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.78077e-05,0.09161,-9.48341e-05,5.51684e-08,-1.3345e-11,34887.8,35.2072], Tmin=(100,'K'), Tmax=(986.245,'K')), NASAPolynomial(coeffs=[12.3957,0.0413358,-1.83714e-05,3.4827e-09,-2.43493e-13,32442.8,-24.4236], Tmin=(986.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])(C=CC=C)C[O](27569)',
    structure = SMILES('[CH2]C([O])(C=CC=C)C[O]'),
    E0 = (322.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.150171,0.0887805,-8.89345e-05,5.06512e-08,-1.21342e-11,38868.5,35.995], Tmin=(100,'K'), Tmax=(989.02,'K')), NASAPolynomial(coeffs=[11.5146,0.0428182,-1.9226e-05,3.66317e-09,-2.56838e-13,36620.6,-18.7066], Tmin=(989.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1([CH][O])C[CH][CH]CO1(27429)',
    structure = SMILES('[CH2]C1([CH][O])C[CH][CH]CO1'),
    E0 = (528.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0712926,0.0876803,-9.53619e-05,6.46106e-08,-1.807e-11,63727.6,32.7515], Tmin=(100,'K'), Tmax=(942.848,'K')), NASAPolynomial(coeffs=[9.64464,0.0422087,-1.5293e-05,2.53202e-09,-1.60892e-13,62138.3,-11.7261], Tmin=(942.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(CCOJ) + radical(RCCJCC) + radical(CCJCO) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1([O])C[CH][CH]CC1(27313)',
    structure = SMILES('[O][CH]C1([O])C[CH][CH]CC1'),
    E0 = (483.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31947,0.0595204,-2.52928e-05,3.08758e-09,7.86062e-14,58156.9,25.5748], Tmin=(100,'K'), Tmax=(2307.33,'K')), NASAPolynomial(coeffs=[33.0775,0.0219868,-1.11564e-05,1.96882e-09,-1.21513e-13,39760.3,-157.639], Tmin=(2307.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclohexane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(cyclohexane) + radical(cyclohexane) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])C[CH][CH]CC1[O](27570)',
    structure = SMILES('[CH2]C1([O])C[CH][CH]CC1[O]'),
    E0 = (506.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09152,0.0690046,-3.81938e-05,9.36977e-09,-8.93154e-13,61070.7,28.2201], Tmin=(100,'K'), Tmax=(2341.08,'K')), NASAPolynomial(coeffs=[22.9265,0.031698,-1.42909e-05,2.56314e-09,-1.66303e-13,50846.9,-95.6945], Tmin=(2341.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclohexane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(cyclohexane) + radical(cyclohexane) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])[CH]OC[CH][CH]C1(27361)',
    structure = SMILES('[CH2]C1([O])[CH]OC[CH][CH]C1'),
    E0 = (558.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.629605,0.0873683,-7.71859e-05,3.75327e-08,-7.19495e-12,67316.4,33.8769], Tmin=(100,'K'), Tmax=(1402.35,'K')), NASAPolynomial(coeffs=[17.1059,0.0297517,-8.03893e-06,1.08661e-09,-6.04502e-14,63033.2,-55.2194], Tmin=(1402.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(oxepane) + radical(CC(C)2OJ) + radical(RCCJCC) + radical(CCJCO) + radical(CCsJOCs) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C([O])([CH][O])C=CC=C(25918)',
    structure = SMILES('[CH2]C([O])([CH][O])C=CC=C'),
    E0 = (502.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,180,180,180,180,888.224,1337.94,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15911,'amu*angstrom^2'), symmetry=1, barrier=(3.65825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15911,'amu*angstrom^2'), symmetry=1, barrier=(3.65825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15911,'amu*angstrom^2'), symmetry=1, barrier=(3.65825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15911,'amu*angstrom^2'), symmetry=1, barrier=(3.65825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (124.137,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.23115,0.0988325,-0.000124718,8.89607e-08,-2.59787e-11,60565.3,36.9131], Tmin=(100,'K'), Tmax=(830.598,'K')), NASAPolynomial(coeffs=[12.0155,0.0398505,-1.81931e-05,3.45371e-09,-2.40173e-13,58531,-19.8959], Tmin=(830.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])([CH][O])[CH]CC=C(27462)',
    structure = SMILES('[CH2]C([O])([CH][O])[CH]CC=C'),
    E0 = (612.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.798106,0.114083,-0.000167284,1.39063e-07,-4.65028e-11,73776,39.6706], Tmin=(100,'K'), Tmax=(810.329,'K')), NASAPolynomial(coeffs=[11.2752,0.0447876,-2.10604e-05,3.99361e-09,-2.74952e-13,72137.8,-14.0722], Tmin=(810.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCJCO) + radical(CCsJOH) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C([O])([CH][O])CC[C]=C(27571)',
    structure = SMILES('[CH2]C([O])([CH][O])CC[C]=C'),
    E0 = (649.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06294,0.124373,-0.000202535,1.80425e-07,-6.2506e-11,78344.4,38.5957], Tmin=(100,'K'), Tmax=(849.103,'K')), NASAPolynomial(coeffs=[9.69868,0.0482686,-2.32063e-05,4.38583e-09,-2.9893e-13,77432.7,-6.169], Tmin=(849.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCCC([CH2])([O])[CH][O](27572)',
    structure = SMILES('[CH]=CCCC([CH2])([O])[CH][O]'),
    E0 = (659.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09584,0.123189,-0.000193957,1.67841e-07,-5.71021e-11,79460.4,38.6442], Tmin=(100,'K'), Tmax=(839.069,'K')), NASAPolynomial(coeffs=[11.1166,0.0459596,-2.19101e-05,4.14201e-09,-2.83081e-13,78080.2,-14.1434], Tmin=(839.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]CC([CH2])(O)[CH][O](27573)',
    structure = SMILES('[CH]C=CCC([CH2])(O)[CH][O]'),
    E0 = (538.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1627,0.123402,-0.000182216,1.53193e-07,-5.13577e-11,64920.7,39.0914], Tmin=(100,'K'), Tmax=(839.188,'K')), NASAPolynomial(coeffs=[10.8409,0.0503121,-2.31984e-05,4.32502e-09,-2.93898e-13,63465,-13.3842], Tmin=(839.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]CC(C)([O])[CH][O](27574)',
    structure = SMILES('[CH]C=CCC(C)([O])[CH][O]'),
    E0 = (555.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.562237,0.10796,-0.000140098,1.12486e-07,-3.78666e-11,66933.9,37.1624], Tmin=(100,'K'), Tmax=(778.391,'K')), NASAPolynomial(coeffs=[8.91763,0.053719,-2.49245e-05,4.72323e-09,-3.26627e-13,65625.5,-5.12224], Tmin=(778.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]CC([CH2])([O])C[O](27575)',
    structure = SMILES('[CH]C=CCC([CH2])([O])C[O]'),
    E0 = (586.932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.636041,0.11021,-0.000147184,1.20933e-07,-4.13206e-11,70750.7,37.715], Tmin=(100,'K'), Tmax=(788.711,'K')), NASAPolynomial(coeffs=[8.74631,0.0544138,-2.54475e-05,4.8303e-09,-3.33925e-13,69526.2,-3.70291], Tmin=(788.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]CC([CH2])([O])[CH]O(27576)',
    structure = SMILES('[CH]C=CCC([CH2])([O])[CH]O'),
    E0 = (541.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (125.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19446,0.120737,-0.000162184,1.22019e-07,-3.719e-11,65311.4,38.4136], Tmin=(100,'K'), Tmax=(800.558,'K')), NASAPolynomial(coeffs=[13.9735,0.0449476,-2.01725e-05,3.75406e-09,-2.56722e-13,62882.9,-31.3887], Tmin=(800.558,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)2OJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(AllylJ2_triplet)"""),
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
    E0 = (548.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (548.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (705.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (718.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1165.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (823.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (939.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1117.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (996.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1040.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (979.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (553.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (550.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (553.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (555.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (556.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (555.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (573.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (573.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (573.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (906.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (773.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (605.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (622.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (610.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (602.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (670.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (617.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (641.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (638.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (796.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (548.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (580.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (548.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (548.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (846.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (941.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (876.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (997.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (997.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (755.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (690.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (667.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (672.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (796.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (627.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (658.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (618.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (706.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (816.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (711.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (601.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (613.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (649.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (634.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (769.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (1122.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (556.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (556.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (556.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (626.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (626.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (626.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (608.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (635.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (614.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (626.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (715.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (673.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (737.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (850.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (816.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (676.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (687.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (738.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (676.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['C=C([O])[CH][O](2850)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=CCC([CH2])([O])C=O(24800)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['C=C[CH]CC[C]([O])[CH][O](25680)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][C]([O])C([O])C[CH]C=C(25655)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2]C=CC[C]([CH2])[CH][O](27537)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]CC([O])=C[O](11358)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][O](1548)', '[CH2][CH][CH]CC(=C)[O](4319)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH]=CCC([CH2])([O])[CH][O](27538)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([O])([CH][O])CC=C[CH2](27539)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2]C=CCC([CH2])([O])[C][O](27540)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]=C[CH]CC([CH2])([O])[CH][O](25909)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=CCC1([CH][O])CO1(27541)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1(C[CH]C=C)OC1[O](24062)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=CCC1([O])CC1[O](27542)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.01367e+13,'s^-1'), n=0.0238333, Ea=(5.50893,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_2H;Ypri_rad_out] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([CH][O])CC=CCO1(25668)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[O][CH]C1([O])CC=CCC1(25682)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_2H] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([O])CC=CCC1[O](25657)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C(O)([CH][O])CC=C=C(25672)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['C=C=CCC(C)([O])[CH][O](25683)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C([O])(C[O])CC=C=C(25662)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C[CH]C([CH2])([O])[CH][O](27543)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C([O])([CH][O])CC1[CH]C1(27544)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1[CH]CC([CH2])([CH][O])O1(27401)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(57.3526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1[CH]CC([O])([CH][O])C1(27235)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.47079e+07,'s^-1'), n=0.909323, Ea=(74.2834,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1[CH]CC([CH2])([O])C1[O](27509)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.5398e+06,'s^-1'), n=1.09629, Ea=(62.0258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1[CH]CC([CH2])([O])[CH]O1(27331)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.18543e+10,'s^-1'), n=0.209443, Ea=(54.3406,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2][CH]C1CC([CH2])([CH][O])O1(25794)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2][CH]C1CC([O])([CH][O])C1(27545)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(69.8718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 65.3 to 69.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2][CH]C1CC([CH2])([O])C1[O](27546)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(93.6028,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 89.1 to 93.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2][CH]C1CC([CH2])([O])[CH]O1(27547)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.34604e+09,'s^-1'), n=0.547362, Ea=(90.0766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R6;doublebond_intra;radadd_intra] for rate rule [R6;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]C([O])([CH][O])CC=C=C(27548)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(T)(20)', 'C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(105.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;Y_1centerbirad] for rate rule [CO_O;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.8 to 105.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH][O](1548)', 'C=C[CH]CC(=C)[O](4317)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.75465,'m^3/(mol*s)'), n=1.70829, Ea=(22.0129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [CO_O;YJ] for rate rule [CO_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(4)', '[CH2]C=CCC(=C)[CH][O](22101)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(55.7026,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 55.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([O])[CH][O](2850)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.000357827,'m^3/(mol*s)'), n=2.74787, Ea=(252.375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-CdHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 251.3 to 252.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C]([O])[CH][O](26240)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.23391e+07,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C[CH2](16804)', '[CH2]C([CH2])([O])[CH][O](26275)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.62e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', '[CH2]C=C[CH]C([CH2])([O])[CH][O](27549)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(3)', '[CH2]C=[C]CC([CH2])([O])[CH][O](27550)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', '[CH2][C]=CCC([CH2])([O])[CH][O](27551)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C([O])([CH][O])CC=[C]C(27552)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=C[CH]C([CH2])(O)[CH][O](27553)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.07519e+07,'s^-1'), n=1.60667, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=C[CH]C(C)([O])[CH][O](27554)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 85 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=C[CH]C([CH2])([O])C[O](27555)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])([CH][O])C[C]=CC(27556)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C=C[CH]C([CH2])([O])[CH]O(27557)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(608922,'s^-1'), n=1.98643, Ea=(79.1566,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4HJ_1;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C=[C]CC([CH2])(O)[CH][O](27558)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C=[C]CC(C)([O])[CH][O](27559)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C=[C]CC([CH2])([O])C[O](27560)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.572e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C([O])([CH][O])[CH]C=CC(27561)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(7.15542e+07,'s^-1'), n=1, Ea=(267.985,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SDS;C_rad_out_single;Cs_H_out_H/NonDeC] + [R4H_SDS;C_rad_out_2H;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C=[C]CC([CH2])([O])[CH]O(27562)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][C]=CCC([CH2])(O)[CH][O](27563)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][C]=CCC(C)([O])[CH][O](27564)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2][C]=CCC([CH2])([O])C[O](27565)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(760143,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;XH_out]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2][C]=CCC([CH2])([O])[CH]O(27566)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_4;Cd_rad_out;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([CH2])([O])[CH][O](24710)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C2H3(60)', '[CH]CC([CH2])([O])[CH][O](26278)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([CH][O])CC(C=C)O1(25666)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['C=CC1CC([O])([CH][O])C1(27229)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([O])CC(C=C)C1[O](27500)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C(O)([CH][O])C=CC=C(27567)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['C=CC=CC(C)([O])[CH][O](27568)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C([O])(C=CC=C)C[O](27569)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([CH][O])C[CH][CH]CO1(27429)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[O][CH]C1([O])C[CH][CH]CC1(27313)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([O])C[CH][CH]CC1[O](27570)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(5.52313e+07,'s^-1'), n=0.767814, Ea=(66.7931,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH2]C1([O])[CH]OC[CH][CH]C1(27361)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(2.76476e+07,'s^-1'), n=0.815689, Ea=(78.927,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction68',
    reactants = ['H(3)', '[CH2]C([O])([CH][O])C=CC=C(25918)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2][C]([O])[CH][O](26240)', 'butadiene13(1350)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]C([O])([CH][O])[CH]CC=C(27462)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH2]C([O])([CH][O])CC[C]=C(27571)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH]=CCCC([CH2])([O])[CH][O](27572)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH]=C[CH]CC([CH2])(O)[CH][O](27573)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH]=C[CH]CC(C)([O])[CH][O](27574)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH]=C[CH]CC([CH2])([O])C[O](27575)'],
    products = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH2]C=CCC([CH2])([O])[CH][O](24788)'],
    products = ['[CH]=C[CH]CC([CH2])([O])[CH]O(27576)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R7Hall;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #6312',
    isomers = [
        '[CH2]C=CCC([CH2])([O])[CH][O](24788)',
    ],
    reactants = [
        ('C=C([O])[CH][O](2850)', 'butadiene13(1350)'),
        ('C=C([O])[CH][O](2850)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #6312',
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


species(
    label = '[CH2]C(C)C([CH2])C=C(3790)',
    structure = SMILES('[CH2]C(C)C([CH2])C=C'),
    E0 = (308.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,637.459,2654.7],'cm^-1')),
        HinderedRotor(inertia=(3.3788,'amu*angstrom^2'), symmetry=1, barrier=(77.6852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247981,'amu*angstrom^2'), symmetry=1, barrier=(14.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00280478,'amu*angstrom^2'), symmetry=1, barrier=(14.0303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611635,'amu*angstrom^2'), symmetry=1, barrier=(14.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.37535,'amu*angstrom^2'), symmetry=1, barrier=(77.606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611916,0.0645992,-2.09321e-05,-1.73983e-08,1.16313e-11,37251.2,32.0101], Tmin=(100,'K'), Tmax=(956.081,'K')), NASAPolynomial(coeffs=[13.0911,0.0357297,-1.2257e-05,2.08623e-09,-1.4024e-13,33798.3,-33.2129], Tmin=(956.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C3H6(27)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,-487.139,-4.54469], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'allyl(82)',
    structure = SMILES('[CH2]C=C'),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29615,0.00579223,4.3392e-05,-5.9989e-08,2.33815e-11,18908.2,9.01994], Tmin=(100,'K'), Tmax=(942.181,'K')), NASAPolynomial(coeffs=[8.06863,0.0101837,-2.84795e-06,5.0088e-10,-3.79628e-14,16914.7,-19.5272], Tmin=(942.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'm1_allyl(186)',
    structure = SMILES('C=C[CH]C'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800698,'amu*angstrom^2'), symmetry=1, barrier=(1.84096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.062381,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]CC([CH2])C=C(3869)',
    structure = SMILES('[CH2]CC([CH2])C=C'),
    E0 = (338.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,3819.65,3819.8],'cm^-1')),
        HinderedRotor(inertia=(0.918967,'amu*angstrom^2'), symmetry=1, barrier=(21.1289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919163,'amu*angstrom^2'), symmetry=1, barrier=(21.1334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.78395,'amu*angstrom^2'), symmetry=1, barrier=(87.0006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127208,'amu*angstrom^2'), symmetry=1, barrier=(86.9922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.55,'J/mol'), sigma=(5.8838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.98 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28275,0.051164,-1.4278e-05,-1.45061e-08,8.61961e-12,40790.4,27.6427], Tmin=(100,'K'), Tmax=(1000.13,'K')), NASAPolynomial(coeffs=[11.0728,0.0305554,-1.11858e-05,1.97483e-09,-1.35046e-13,37904.6,-24.2271], Tmin=(1000.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C=C)C[CH]C(3951)',
    structure = SMILES('[CH2]C(C=C)C[CH]C'),
    E0 = (303.694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69878,0.0664879,-4.02964e-05,1.30928e-08,-1.79722e-12,36649.6,32.9385], Tmin=(100,'K'), Tmax=(1622.71,'K')), NASAPolynomial(coeffs=[12.0454,0.0385183,-1.44419e-05,2.47082e-09,-1.60772e-13,32967.1,-27.2955], Tmin=(1622.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]CC)C=C(4977)',
    structure = SMILES('[CH2]C([CH]CC)C=C'),
    E0 = (303.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558145,0.0665336,-3.70022e-05,8.38567e-09,-1.53865e-13,36668.9,33.2709], Tmin=(100,'K'), Tmax=(1271.49,'K')), NASAPolynomial(coeffs=[12.6063,0.0381477,-1.47419e-05,2.60079e-09,-1.73869e-13,32835.8,-30.7736], Tmin=(1271.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C[CH]C=C(3936)',
    structure = SMILES('[CH2]C=CCC([CH2])C'),
    E0 = (247.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0795356,'amu*angstrom^2'), symmetry=1, barrier=(1.82868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00784389,'amu*angstrom^2'), symmetry=1, barrier=(1.95248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.082228,'amu*angstrom^2'), symmetry=1, barrier=(1.89058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0747745,'amu*angstrom^2'), symmetry=1, barrier=(18.8182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.42126,'amu*angstrom^2'), symmetry=1, barrier=(78.6616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.572947,0.064088,-1.83159e-05,-1.83223e-08,1.09052e-11,29892.5,30.4562], Tmin=(100,'K'), Tmax=(1006.91,'K')), NASAPolynomial(coeffs=[13.5259,0.0368006,-1.36698e-05,2.43925e-09,-1.68027e-13,26058.8,-38.2074], Tmin=(1006.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)[CH]CC=C(4978)',
    structure = SMILES('[CH2]C(C)[CH]CC=C'),
    E0 = (306.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.454728,0.0678108,-4.10446e-05,1.27254e-08,-1.6161e-12,36947.6,33.6235], Tmin=(100,'K'), Tmax=(1791.32,'K')), NASAPolynomial(coeffs=[15.6952,0.0337788,-1.2547e-05,2.11946e-09,-1.35909e-13,31487.5,-48.7871], Tmin=(1791.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH]C)C=C(3803)',
    structure = SMILES('[CH2]C([CH]C)C=C'),
    E0 = (327.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1026.84,3543.87],'cm^-1')),
        HinderedRotor(inertia=(0.235199,'amu*angstrom^2'), symmetry=1, barrier=(5.40768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235218,'amu*angstrom^2'), symmetry=1, barrier=(5.40812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271982,'amu*angstrom^2'), symmetry=1, barrier=(20.3613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.76263,'amu*angstrom^2'), symmetry=1, barrier=(86.5103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39619,0.0496388,-1.64699e-05,-7.33995e-09,4.82136e-12,39497.8,28.008], Tmin=(100,'K'), Tmax=(1082.86,'K')), NASAPolynomial(coeffs=[9.7195,0.0328979,-1.26799e-05,2.27041e-09,-1.54839e-13,36874.1,-16.601], Tmin=(1082.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=CC([CH2])C(532)',
    structure = SMILES('[CH2]C=CC([CH2])C'),
    E0 = (272.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.00511567,'amu*angstrom^2'), symmetry=1, barrier=(4.60508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679991,'amu*angstrom^2'), symmetry=1, barrier=(15.6343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0173188,'amu*angstrom^2'), symmetry=1, barrier=(15.6385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.30708,'amu*angstrom^2'), symmetry=1, barrier=(76.0364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32269,0.0477996,1.34752e-06,-3.4245e-08,1.64564e-11,32854.3,25.5549], Tmin=(100,'K'), Tmax=(967.436,'K')), NASAPolynomial(coeffs=[12.2804,0.0288314,-1.00791e-05,1.77001e-09,-1.22476e-13,29501.6,-33.3169], Tmin=(967.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = '[CH]C(C)C([CH2])C=C(4979)',
    structure = SMILES('[CH]C(C)C([CH2])C=C'),
    E0 = (551.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51754,0.0662538,-2.89555e-05,-8.87045e-09,8.26124e-12,66496.9,31.5712], Tmin=(100,'K'), Tmax=(992.953,'K')), NASAPolynomial(coeffs=[14.293,0.0328624,-1.19004e-05,2.09517e-09,-1.43482e-13,62671.7,-40.2765], Tmin=(992.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C)C([CH2])C(4980)',
    structure = SMILES('[CH]C(C=C)C([CH2])C'),
    E0 = (551.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51754,0.0662538,-2.89555e-05,-8.87045e-09,8.26124e-12,66496.9,31.5712], Tmin=(100,'K'), Tmax=(992.953,'K')), NASAPolynomial(coeffs=[14.293,0.0328624,-1.19004e-05,2.09517e-09,-1.43482e-13,62671.7,-40.2765], Tmin=(992.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=CC1CCC1C(4981)',
    structure = SMILES('C=CC1CCC1C'),
    E0 = (53.3011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50542,0.0308318,8.53225e-05,-1.29014e-07,5.04011e-11,6522.1,24.1147], Tmin=(100,'K'), Tmax=(966.994,'K')), NASAPolynomial(coeffs=[15.5153,0.033087,-1.15692e-05,2.17273e-09,-1.61297e-13,997.679,-57.5596], Tmin=(966.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.3011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=CC(C)C(=C)C(4982)',
    structure = SMILES('C=CC(C)C(=C)C'),
    E0 = (16.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607449,0.070627,-4.57496e-05,1.53761e-08,-2.1561e-12,2064.29,25.3447], Tmin=(100,'K'), Tmax=(1597.76,'K')), NASAPolynomial(coeffs=[13.4412,0.0384975,-1.55859e-05,2.79019e-09,-1.86787e-13,-2036.73,-42.5845], Tmin=(1597.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=C)C(C)C(4983)',
    structure = SMILES('C=CC(=C)C(C)C'),
    E0 = (1.93299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.490126,0.0615741,-7.04723e-07,-4.38923e-08,2.15655e-11,372.822,25.4118], Tmin=(100,'K'), Tmax=(975.04,'K')), NASAPolynomial(coeffs=[16.8015,0.0315102,-1.11471e-05,2.0101e-09,-1.42514e-13,-4559.78,-61.8518], Tmin=(975.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.93299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C[C]([CH2])C([CH2])C(4984)',
    structure = SMILES('[CH2]C[C]([CH2])C([CH2])C'),
    E0 = (577.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,365.412,3970.37],'cm^-1')),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.584348,0.0789203,-8.25247e-05,5.9668e-08,-1.8791e-11,69534.2,34.5259], Tmin=(100,'K'), Tmax=(840.351,'K')), NASAPolynomial(coeffs=[6.16113,0.0479985,-1.95178e-05,3.48565e-09,-2.33261e-13,68751.4,9.51069], Tmin=(840.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])[C]([CH2])C(4985)',
    structure = SMILES('[CH2]CC([CH2])[C]([CH2])C'),
    E0 = (577.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,365.412,3970.37],'cm^-1')),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0034047,'amu*angstrom^2'), symmetry=1, barrier=(0.20097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.584348,0.0789203,-8.25247e-05,5.9668e-08,-1.8791e-11,69534.2,34.5259], Tmin=(100,'K'), Tmax=(840.351,'K')), NASAPolynomial(coeffs=[6.16113,0.0479985,-1.95178e-05,3.48565e-09,-2.33261e-13,68751.4,9.51069], Tmin=(840.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C([CH2])[CH]C(3924)',
    structure = SMILES('[CH2][C](C)C([CH2])[CH]C'),
    E0 = (566.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,917.845,3925.58],'cm^-1')),
        HinderedRotor(inertia=(0.0159029,'amu*angstrom^2'), symmetry=1, barrier=(7.56243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159029,'amu*angstrom^2'), symmetry=1, barrier=(7.56243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159029,'amu*angstrom^2'), symmetry=1, barrier=(7.56243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159029,'amu*angstrom^2'), symmetry=1, barrier=(7.56243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159029,'amu*angstrom^2'), symmetry=1, barrier=(7.56243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0159029,'amu*angstrom^2'), symmetry=1, barrier=(7.56243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.743385,0.076926,-8.34602e-05,6.58885e-08,-2.25483e-11,68239.5,34.7236], Tmin=(100,'K'), Tmax=(838.257,'K')), NASAPolynomial(coeffs=[4.24393,0.0512716,-2.15373e-05,3.90338e-09,-2.63066e-13,67967.1,20.3288], Tmin=(838.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])[CH2](571)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])[CH2]'),
    E0 = (596.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,268.107,2267.26],'cm^-1')),
        HinderedRotor(inertia=(0.0564368,'amu*angstrom^2'), symmetry=1, barrier=(3.0641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0564368,'amu*angstrom^2'), symmetry=1, barrier=(3.0641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0564368,'amu*angstrom^2'), symmetry=1, barrier=(3.0641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0564368,'amu*angstrom^2'), symmetry=1, barrier=(3.0641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0564368,'amu*angstrom^2'), symmetry=1, barrier=(3.0641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0564368,'amu*angstrom^2'), symmetry=1, barrier=(3.0641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140897,0.0776731,-6.40399e-05,3.08497e-08,-6.14544e-12,71924.3,36.0654], Tmin=(100,'K'), Tmax=(1201.59,'K')), NASAPolynomial(coeffs=[12.7062,0.0358433,-1.1821e-05,1.87713e-09,-1.17384e-13,68904.7,-26.8625], Tmin=(1201.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])[CH]C(570)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[CH]C'),
    E0 = (586.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,304.417,3135.38],'cm^-1')),
        HinderedRotor(inertia=(0.0305351,'amu*angstrom^2'), symmetry=1, barrier=(3.73581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0305351,'amu*angstrom^2'), symmetry=1, barrier=(3.73581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0305351,'amu*angstrom^2'), symmetry=1, barrier=(3.73581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0305351,'amu*angstrom^2'), symmetry=1, barrier=(3.73581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0305351,'amu*angstrom^2'), symmetry=1, barrier=(3.73581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0305351,'amu*angstrom^2'), symmetry=1, barrier=(3.73581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3395.46,'J/mol'), sigma=(6.43099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.36 K, Pc=28.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478804,0.0735456,-5.73644e-05,2.68911e-08,-5.3752e-12,70621.9,35.6228], Tmin=(100,'K'), Tmax=(1173.54,'K')), NASAPolynomial(coeffs=[10.3807,0.0397952,-1.42252e-05,2.38451e-09,-1.54548e-13,68297.9,-13.7325], Tmin=(1173.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C1[CH]CC1(4986)',
    structure = SMILES('[CH2]C(C)C1[CH]CC1'),
    E0 = (317.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37408,0.0412413,4.3258e-05,-7.81611e-08,3.09923e-11,38280,28.8997], Tmin=(100,'K'), Tmax=(982.333,'K')), NASAPolynomial(coeffs=[11.9123,0.0386111,-1.42337e-05,2.59878e-09,-1.84108e-13,34266.1,-31.6458], Tmin=(982.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]CCC1C(4987)',
    structure = SMILES('[CH2]C1[CH]CCC1C'),
    E0 = (242.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76703,0.0279975,8.33815e-05,-1.22355e-07,4.75591e-11,29242.3,27.1969], Tmin=(100,'K'), Tmax=(956.507,'K')), NASAPolynomial(coeffs=[12.6288,0.0356856,-1.19635e-05,2.14885e-09,-1.54716e-13,24734.9,-37.4218], Tmin=(956.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C1CC1[CH2](4988)',
    structure = SMILES('[CH2]C(C)C1CC1[CH2]'),
    E0 = (329.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00473,0.0482094,3.51012e-05,-8.17472e-08,3.60841e-11,39812.9,28.6728], Tmin=(100,'K'), Tmax=(932.012,'K')), NASAPolynomial(coeffs=[15.5381,0.0307251,-9.00599e-06,1.48059e-09,-1.02573e-13,35154.2,-50.8786], Tmin=(932.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC(C)C1[CH2](4101)',
    structure = SMILES('[CH2]C1CC(C)C1[CH2]'),
    E0 = (325.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29058,0.0392127,6.12846e-05,-1.07404e-07,4.46018e-11,39273.5,26.9676], Tmin=(100,'K'), Tmax=(937.538,'K')), NASAPolynomial(coeffs=[15.1048,0.0318576,-9.47711e-06,1.59846e-09,-1.13149e-13,34416.3,-50.8771], Tmin=(937.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C(=C)C(4989)',
    structure = SMILES('[CH2]C(C=C)C(=C)C'),
    E0 = (221.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,210.202,210.202],'cm^-1')),
        HinderedRotor(inertia=(0.218669,'amu*angstrom^2'), symmetry=1, barrier=(6.8563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00381517,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218673,'amu*angstrom^2'), symmetry=1, barrier=(6.8564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.669382,'amu*angstrom^2'), symmetry=1, barrier=(20.9871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726623,0.0711765,-5.47003e-05,2.40787e-08,-4.54742e-12,26723.5,26.6972], Tmin=(100,'K'), Tmax=(1217.8,'K')), NASAPolynomial(coeffs=[10.2983,0.0397368,-1.59746e-05,2.87863e-09,-1.95243e-13,24392.3,-21.3668], Tmin=(1217.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C(=C)C=C(4990)',
    structure = SMILES('[CH2]C(C)C(=C)C=C'),
    E0 = (207.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,385.721],'cm^-1')),
        HinderedRotor(inertia=(0.690583,'amu*angstrom^2'), symmetry=1, barrier=(15.8779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690591,'amu*angstrom^2'), symmetry=1, barrier=(15.878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0377624,'amu*angstrom^2'), symmetry=1, barrier=(15.8781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037761,'amu*angstrom^2'), symmetry=1, barrier=(15.8778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526603,0.063039,-1.2415e-05,-3.24748e-08,1.84538e-11,25035.4,27.7558], Tmin=(100,'K'), Tmax=(947.869,'K')), NASAPolynomial(coeffs=[16.4982,0.0282957,-9.11282e-06,1.55024e-09,-1.06947e-13,20540.6,-56.1814], Tmin=(947.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
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
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([532.913,1391.12,1391.12,2779.21,3448.45,3448.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.948,0.0008276,8.34932e-06,-9.82634e-09,3.80104e-12,16425.4,0.336655], Tmin=(100,'K'), Tmax=(660.467,'K')), NASAPolynomial(coeffs=[3.2217,0.00522646,-1.64125e-06,2.58225e-10,-1.62579e-14,16521.3,3.53938], Tmin=(660.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C=C)C=C(4888)',
    structure = SMILES('[CH2]C(C=C)C=C'),
    E0 = (260.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.164237,'amu*angstrom^2'), symmetry=1, barrier=(3.77612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16638,'amu*angstrom^2'), symmetry=1, barrier=(3.8254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14555,'amu*angstrom^2'), symmetry=1, barrier=(26.3384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22057,0.0552061,-3.78148e-05,1.41306e-08,-2.20454e-12,31407.2,23.4798], Tmin=(100,'K'), Tmax=(1478.15,'K')), NASAPolynomial(coeffs=[11.1866,0.0282371,-1.0447e-05,1.78734e-09,-1.1691e-13,28461,-28.495], Tmin=(1478.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = 'C3H6(T)(28)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0090964,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93779,0.019099,4.26879e-06,-1.44878e-08,5.7496e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93905,0.0171893,-6.69156e-06,1.21547e-09,-8.39803e-14,33151.2,-4.14862], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(C)C=C(109)',
    structure = SMILES('[CH2]C(C)C=C'),
    E0 = (156.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,663.668],'cm^-1')),
        HinderedRotor(inertia=(0.00204431,'amu*angstrom^2'), symmetry=1, barrier=(15.8768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690941,'amu*angstrom^2'), symmetry=1, barrier=(15.8861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.57515,'amu*angstrom^2'), symmetry=1, barrier=(82.1997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99035,0.0344738,1.41962e-05,-4.03588e-08,1.76162e-11,18940.8,21.072], Tmin=(100,'K'), Tmax=(956.909,'K')), NASAPolynomial(coeffs=[9.88824,0.0252392,-8.60353e-06,1.49487e-09,-1.03137e-13,16340.6,-22.3716], Tmin=(956.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])C=C(4950)',
    structure = SMILES('[CH2][CH]C([CH2])C=C'),
    E0 = (532.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,825.24,825.296],'cm^-1')),
        HinderedRotor(inertia=(0.112603,'amu*angstrom^2'), symmetry=1, barrier=(2.58896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112591,'amu*angstrom^2'), symmetry=1, barrier=(2.58869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0053604,'amu*angstrom^2'), symmetry=1, barrier=(2.59117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112621,'amu*angstrom^2'), symmetry=1, barrier=(2.58939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36168,0.0511269,-2.91252e-05,7.14084e-09,-2.95895e-13,64183.7,29.9243], Tmin=(100,'K'), Tmax=(1256.84,'K')), NASAPolynomial(coeffs=[10.4271,0.0291586,-1.11215e-05,1.94845e-09,-1.29797e-13,61361.3,-18.046], Tmin=(1256.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C)C([CH2])C=C(4991)',
    structure = SMILES('[CH2][C](C)C([CH2])C=C'),
    E0 = (494.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,1562.96],'cm^-1')),
        HinderedRotor(inertia=(0.0899948,'amu*angstrom^2'), symmetry=1, barrier=(2.06916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0899913,'amu*angstrom^2'), symmetry=1, barrier=(2.06908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866916,'amu*angstrom^2'), symmetry=1, barrier=(1.99321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889219,'amu*angstrom^2'), symmetry=1, barrier=(2.04449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0862828,'amu*angstrom^2'), symmetry=1, barrier=(1.98381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907313,0.0694603,-5.84523e-05,3.19267e-08,-7.76034e-12,59531.6,32.0914], Tmin=(100,'K'), Tmax=(956.406,'K')), NASAPolynomial(coeffs=[7.23037,0.0430144,-1.6974e-05,3.01321e-09,-2.02264e-13,58322.2,1.86836], Tmin=(956.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(111)',
    structure = SMILES('[CH2][CH]C([CH2])C'),
    E0 = (431.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1088.38],'cm^-1')),
        HinderedRotor(inertia=(0.00507348,'amu*angstrom^2'), symmetry=1, barrier=(4.26448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00507018,'amu*angstrom^2'), symmetry=1, barrier=(4.2639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185504,'amu*angstrom^2'), symmetry=1, barrier=(4.2651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185366,'amu*angstrom^2'), symmetry=1, barrier=(4.26194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72887,0.0438928,-2.37827e-05,6.55025e-09,-7.41727e-13,51935.3,25.8655], Tmin=(100,'K'), Tmax=(1961.54,'K')), NASAPolynomial(coeffs=[11.3781,0.0242161,-8.73612e-06,1.43645e-09,-8.99776e-14,48149.8,-27.1877], Tmin=(1961.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])C(4992)',
    structure = SMILES('[CH2]C=C([CH2])C([CH2])C'),
    E0 = (388.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.702586,'amu*angstrom^2'), symmetry=1, barrier=(16.1538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187656,'amu*angstrom^2'), symmetry=1, barrier=(4.31459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0219808,'amu*angstrom^2'), symmetry=1, barrier=(16.1538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111642,'amu*angstrom^2'), symmetry=1, barrier=(82.0812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111592,'amu*angstrom^2'), symmetry=1, barrier=(82.07,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486903,0.065208,-2.02721e-05,-2.10962e-08,1.33404e-11,46912.4,29.592], Tmin=(100,'K'), Tmax=(971.957,'K')), NASAPolynomial(coeffs=[15.41,0.0315193,-1.10698e-05,1.94083e-09,-1.33894e-13,42701.8,-48.716], Tmin=(971.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])C=C(4122)',
    structure = SMILES('[CH2]C([CH2])C([CH2])C=C'),
    E0 = (513.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,405.633,1731.94],'cm^-1')),
        HinderedRotor(inertia=(0.601374,'amu*angstrom^2'), symmetry=1, barrier=(70.2164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0863067,'amu*angstrom^2'), symmetry=1, barrier=(10.0769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00473428,'amu*angstrom^2'), symmetry=1, barrier=(10.0771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00473413,'amu*angstrom^2'), symmetry=1, barrier=(10.0769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.601372,'amu*angstrom^2'), symmetry=1, barrier=(70.2165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683737,0.0655995,-3.07701e-05,-8.81546e-09,9.93426e-12,61912.3,32.844], Tmin=(100,'K'), Tmax=(899.328,'K')), NASAPolynomial(coeffs=[12.7877,0.0325249,-1.02323e-05,1.62934e-09,-1.04962e-13,58895.6,-28.9343], Tmin=(899.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C([CH2])[C]=C(4993)',
    structure = SMILES('[CH2]C(C)C([CH2])[C]=C'),
    E0 = (546.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.148322,'amu*angstrom^2'), symmetry=1, barrier=(83.1789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0776037,'amu*angstrom^2'), symmetry=1, barrier=(12.0895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525751,'amu*angstrom^2'), symmetry=1, barrier=(12.088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0097807,'amu*angstrom^2'), symmetry=1, barrier=(83.1847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00978098,'amu*angstrom^2'), symmetry=1, barrier=(83.1874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.327805,0.0720972,-5.55306e-05,2.40786e-08,-4.30472e-12,65865.5,33.5069], Tmin=(100,'K'), Tmax=(1328.31,'K')), NASAPolynomial(coeffs=[13.3671,0.032831,-1.11888e-05,1.82369e-09,-1.16118e-13,62401.5,-33.1021], Tmin=(1328.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C([CH2])C(4994)',
    structure = SMILES('[CH]=CC([CH2])C([CH2])C'),
    E0 = (555.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,375.124],'cm^-1')),
        HinderedRotor(inertia=(0.00119798,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096574,'amu*angstrom^2'), symmetry=1, barrier=(9.64357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788099,'amu*angstrom^2'), symmetry=1, barrier=(78.6975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965744,'amu*angstrom^2'), symmetry=1, barrier=(9.64356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0116482,'amu*angstrom^2'), symmetry=1, barrier=(78.6975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542517,0.0680693,-3.74016e-05,-2.63568e-10,5.82203e-12,66970.8,32.6624], Tmin=(100,'K'), Tmax=(951.595,'K')), NASAPolynomial(coeffs=[13.3172,0.0329208,-1.12368e-05,1.89074e-09,-1.25644e-13,63699.7,-32.7473], Tmin=(951.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)[C](C)C(4995)',
    structure = SMILES('[CH2]C(C=C)[C](C)C'),
    E0 = (288.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00662,0.0665492,-4.22219e-05,1.50729e-08,-2.39037e-12,34862.6,29.2476], Tmin=(100,'K'), Tmax=(1360.7,'K')), NASAPolynomial(coeffs=[8.4455,0.0446811,-1.81147e-05,3.26159e-09,-2.20256e-13,32838.2,-8.93175], Tmin=(1360.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C(C)C([CH2])C(3889)',
    structure = SMILES('[CH2]C(C)[C](C)C=C'),
    E0 = (235.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2458.99,4000],'cm^-1')),
        HinderedRotor(inertia=(1.6275,'amu*angstrom^2'), symmetry=1, barrier=(96.6712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332426,'amu*angstrom^2'), symmetry=1, barrier=(19.7579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00201695,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62593,'amu*angstrom^2'), symmetry=1, barrier=(96.5782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62086,'amu*angstrom^2'), symmetry=1, barrier=(96.6977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85271,0.0577382,-3.72774e-06,-3.16863e-08,1.54826e-11,28451.7,27.6175], Tmin=(100,'K'), Tmax=(981.607,'K')), NASAPolynomial(coeffs=[12.195,0.0378732,-1.36439e-05,2.39939e-09,-1.64408e-13,24955.3,-33.3592], Tmin=(981.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C=C)C(C)C(4996)',
    structure = SMILES('[CH2]C=C([CH2])C(C)C'),
    E0 = (183.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451556,0.0637135,-8.37192e-06,-3.29051e-08,1.66924e-11,22249.8,27.2449], Tmin=(100,'K'), Tmax=(1000.42,'K')), NASAPolynomial(coeffs=[15.7826,0.0346213,-1.30414e-05,2.38626e-09,-1.68288e-13,17570.6,-54.78], Tmin=(1000.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)C(=C)C(3890)',
    structure = SMILES('[CH2][CH]C(C)C(=C)C'),
    E0 = (292.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,507.101,1106.29],'cm^-1')),
        HinderedRotor(inertia=(0.185042,'amu*angstrom^2'), symmetry=1, barrier=(4.62729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00523978,'amu*angstrom^2'), symmetry=1, barrier=(4.50474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121115,'amu*angstrom^2'), symmetry=1, barrier=(4.25726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111712,'amu*angstrom^2'), symmetry=1, barrier=(4.38171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348649,'amu*angstrom^2'), symmetry=1, barrier=(14.6451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622805,0.0678364,-4.07379e-05,1.2129e-08,-1.47071e-12,35349.9,31.5611], Tmin=(100,'K'), Tmax=(1850.54,'K')), NASAPolynomial(coeffs=[16.0213,0.0345519,-1.37583e-05,2.40943e-09,-1.57642e-13,29650.8,-52.2053], Tmin=(1850.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)C(C)[C]=C(4997)',
    structure = SMILES('[CH2]C(C)C(C)[C]=C'),
    E0 = (341.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,366.534,366.569],'cm^-1')),
        HinderedRotor(inertia=(0.0610606,'amu*angstrom^2'), symmetry=1, barrier=(5.8215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00125478,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610656,'amu*angstrom^2'), symmetry=1, barrier=(5.82148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610579,'amu*angstrom^2'), symmetry=1, barrier=(5.82134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256184,'amu*angstrom^2'), symmetry=1, barrier=(24.4235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479795,0.0684426,-3.63175e-05,3.14762e-09,2.78095e-12,41194.8,31.1778], Tmin=(100,'K'), Tmax=(1065.71,'K')), NASAPolynomial(coeffs=[12.6903,0.0376396,-1.41132e-05,2.48888e-09,-1.68415e-13,37738.8,-32.5115], Tmin=(1065.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])C(C)C=C(1262)',
    structure = SMILES('[CH2]C([CH2])C(C)C=C'),
    E0 = (308.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,637.459,2654.7],'cm^-1')),
        HinderedRotor(inertia=(3.3788,'amu*angstrom^2'), symmetry=1, barrier=(77.6852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247981,'amu*angstrom^2'), symmetry=1, barrier=(14.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00280478,'amu*angstrom^2'), symmetry=1, barrier=(14.0303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611635,'amu*angstrom^2'), symmetry=1, barrier=(14.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.37535,'amu*angstrom^2'), symmetry=1, barrier=(77.606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611916,0.0645992,-2.09321e-05,-1.73983e-08,1.16313e-11,37251.2,31.3169], Tmin=(100,'K'), Tmax=(956.081,'K')), NASAPolynomial(coeffs=[13.0911,0.0357297,-1.2257e-05,2.08623e-09,-1.4024e-13,33798.3,-33.906], Tmin=(956.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)C(4998)',
    structure = SMILES('[CH2]C([C]=C)C(C)C'),
    E0 = (341.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,366.534,366.569],'cm^-1')),
        HinderedRotor(inertia=(0.0610606,'amu*angstrom^2'), symmetry=1, barrier=(5.8215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00125478,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610656,'amu*angstrom^2'), symmetry=1, barrier=(5.82148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610579,'amu*angstrom^2'), symmetry=1, barrier=(5.82134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256184,'amu*angstrom^2'), symmetry=1, barrier=(24.4235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479795,0.0684426,-3.63175e-05,3.14762e-09,2.78095e-12,41194.8,30.4847], Tmin=(100,'K'), Tmax=(1065.71,'K')), NASAPolynomial(coeffs=[12.6903,0.0376396,-1.41132e-05,2.48888e-09,-1.68415e-13,37738.8,-33.2046], Tmin=(1065.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C([CH2])C(4999)',
    structure = SMILES('[CH]=CC(C)C([CH2])C'),
    E0 = (350.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,470.943],'cm^-1')),
        HinderedRotor(inertia=(0.0490184,'amu*angstrom^2'), symmetry=1, barrier=(7.71484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0490168,'amu*angstrom^2'), symmetry=1, barrier=(7.715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0490208,'amu*angstrom^2'), symmetry=1, barrier=(7.71488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.049015,'amu*angstrom^2'), symmetry=1, barrier=(7.71491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172651,'amu*angstrom^2'), symmetry=1, barrier=(27.1711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483929,0.0668613,-2.65614e-05,-1.06134e-08,8.51367e-12,42309.2,31.091], Tmin=(100,'K'), Tmax=(1002.84,'K')), NASAPolynomial(coeffs=[13.7533,0.0359125,-1.31437e-05,2.32069e-09,-1.58743e-13,38542.6,-38.4735], Tmin=(1002.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(C)C(5000)',
    structure = SMILES('[CH]=CC([CH2])C(C)C'),
    E0 = (350.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,470.943],'cm^-1')),
        HinderedRotor(inertia=(0.0490184,'amu*angstrom^2'), symmetry=1, barrier=(7.71484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0490168,'amu*angstrom^2'), symmetry=1, barrier=(7.715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0490208,'amu*angstrom^2'), symmetry=1, barrier=(7.71488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.049015,'amu*angstrom^2'), symmetry=1, barrier=(7.71491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172651,'amu*angstrom^2'), symmetry=1, barrier=(27.1711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483929,0.0668613,-2.65614e-05,-1.06134e-08,8.51367e-12,42309.2,30.3979], Tmin=(100,'K'), Tmax=(1002.84,'K')), NASAPolynomial(coeffs=[13.7533,0.0359125,-1.31437e-05,2.32069e-09,-1.58743e-13,38542.6,-39.1667], Tmin=(1002.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
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
    E0 = (308.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (757.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (468.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (503.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (468.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (468.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (765.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (710.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (763.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (763.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (316.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (372.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (372.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (599.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (640.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (574.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (621.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (611.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (432.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (423.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (345.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (367.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (432.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (418.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (308.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (426.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (403.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (457.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (559.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (669.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (705.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (718.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (600.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (725.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (758.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (767.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (496.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (411.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (438.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (426.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (426.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (493.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (391.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (385.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (394.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (383.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['C3H6(27)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]CC([CH2])C=C(3869)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C(C=C)C[CH]C(3951)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C([CH]CC)C=C(4977)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C(C)C[CH]C=C(3936)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C(C)[CH]CC=C(4978)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2]C([CH]C)C=C(3803)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH2]C=CC([CH2])C(532)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C(C)C([CH2])C=C(4979)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C(C=C)C([CH2])C(4980)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['C=CC1CCC1C(4981)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['C=CC(C)C(=C)C(4982)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['C=CC(=C)C(C)C(4983)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C[C]([CH2])C([CH2])C(4984)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC([CH2])[C]([CH2])C(4985)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C](C)C([CH2])[CH]C(3924)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC([CH2])C([CH2])[CH2](571)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])C([CH2])[CH]C(570)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.27566e+10,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C(C)C1[CH]CC1(4986)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C1[CH]CCC1C(4987)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.12e+09,'s^-1'), n=0.63, Ea=(114.642,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C(C)C1CC1[CH2](4988)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C1CC(C)C1[CH2](4101)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.07e+06,'s^-1'), n=1.46, Ea=(58.9944,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C(C=C)C(=C)C(4989)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C(C)C(=C)C=C(4990)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H6(27)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00337229,'m^3/(mol*s)'), n=2.52599, Ea=(27.9458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 23.4 to 27.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH3(17)', '[CH2]C(C=C)C=C(4888)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(20000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C3H6(T)(28)', 'butadiene13(1350)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C2H3(60)', '[CH2]C(C)C=C(109)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 432 used for Cds-CsH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CsH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C3H6(T)(28)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH3(17)', '[CH2][CH]C([CH2])C=C(4950)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.49263e+08,'m^3/(mol*s)'), n=-0.519467, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00186367249776, var=0.658318197044, Tref=1000.0, N=7, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O
    Total Standard Deviation in ln(k): 1.6312606891
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][C](C)C([CH2])C=C(4991)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C2H3(60)', '[CH2][CH]C([CH2])C(111)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.4834e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2][C](C=C)C([CH2])C(4992)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2]C([CH2])C([CH2])C=C(4122)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.835e+08,'m^3/(mol*s)'), n=0.115342, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH2]C(C)C([CH2])[C]=C(4993)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH]=CC([CH2])C([CH2])C(4994)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['allyl(82)', 'm1_allyl(186)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.71146e+34,'s^-1'), n=-6.74695, Ea=(188.22,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.34291747530997396, var=16.416542173868212, Tref=1000.0, N=7, correlation='Root_1R!H->C_2R!H->C_N-5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_N-5R!H->O
    Total Standard Deviation in ln(k): 8.984253428860972
Exact match found for rate rule [Root_1R!H->C_2R!H->C_N-5R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: Retroene"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2]C(C=C)[C](C)C(4995)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C=C(C)C([CH2])C(3889)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2][C](C=C)C(C)C(4996)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(C)C([CH2])C=C(3790)'],
    products = ['[CH2][CH]C(C)C(=C)C(3890)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(C)C(C)[C]=C(4997)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([CH2])C(C)C=C(1262)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(228000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([C]=C)C(C)C(4998)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(222600,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=CC(C)C([CH2])C(4999)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=CC([CH2])C(C)C(5000)'],
    products = ['[CH2]C(C)C([CH2])C=C(3790)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(816000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1287',
    isomers = [
        '[CH2]C(C)C([CH2])C=C(3790)',
    ],
    reactants = [
        ('C3H6(27)', 'butadiene13(1350)'),
        ('allyl(82)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1287',
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


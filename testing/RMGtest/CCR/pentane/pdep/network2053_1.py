species(
    label = '[CH2][CH]C([CH2])C([CH2])[CH2](4093)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[CH2]'),
    E0 = (791.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,935.446,4000],'cm^-1')),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594522,0.0733338,-6.44468e-05,3.47395e-08,-7.93788e-12,95301.3,36.9958], Tmin=(100,'K'), Tmax=(1040.15,'K')), NASAPolynomial(coeffs=[9.70761,0.0382888,-1.39089e-05,2.34833e-09,-1.52743e-13,93405.4,-7.32843], Tmin=(1040.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(791.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3397.52,'J/mol'), sigma=(6.23225,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.69 K, Pc=31.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683737,0.0655995,-3.07701e-05,-8.81546e-09,9.93426e-12,61912.3,32.844], Tmin=(100,'K'), Tmax=(899.328,'K')), NASAPolynomial(coeffs=[12.7877,0.0325249,-1.02323e-05,1.62934e-09,-1.04962e-13,58895.6,-28.9343], Tmin=(899.328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[CH2](4136)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[CH2]'),
    E0 = (782.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1213.89,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3372.64,'J/mol'), sigma=(6.39812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.80 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898412,0.0752073,-9.10451e-05,7.75073e-08,-2.6973e-11,94275.3,37.1477], Tmin=(100,'K'), Tmax=(883.59,'K')), NASAPolynomial(coeffs=[2.98448,0.0496822,-2.0413e-05,3.61754e-09,-2.38926e-13,94534.4,30.8942], Tmin=(883.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])[CH][CH2](4108)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH][CH2]'),
    E0 = (786.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,555.363,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0136821,'amu*angstrom^2'), symmetry=1, barrier=(14.3727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136821,'amu*angstrom^2'), symmetry=1, barrier=(14.3727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136821,'amu*angstrom^2'), symmetry=1, barrier=(14.3727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136821,'amu*angstrom^2'), symmetry=1, barrier=(14.3727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136821,'amu*angstrom^2'), symmetry=1, barrier=(14.3727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136821,'amu*angstrom^2'), symmetry=1, barrier=(14.3727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.967175,0.0712027,-7.16539e-05,5.33178e-08,-1.80799e-11,94707.7,37.5949], Tmin=(100,'K'), Tmax=(781.218,'K')), NASAPolynomial(coeffs=[4.62312,0.0486708,-2.10702e-05,3.90409e-09,-2.67708e-13,94252.8,21.6044], Tmin=(781.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(786.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH2](502)',
    structure = SMILES('[CH][CH2]'),
    E0 = (557.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1433.18,1433.5],'cm^-1')),
        HinderedRotor(inertia=(0.00559429,'amu*angstrom^2'), symmetry=1, barrier=(8.15686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7298,0.00311091,1.59755e-05,-2.26123e-08,9.32732e-12,67119.1,8.34543], Tmin=(100,'K'), Tmax=(875.199,'K')), NASAPolynomial(coeffs=[4.97018,0.00511238,-6.01253e-07,2.8794e-11,-6.01916e-16,66608.3,0.848354], Tmin=(875.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH2](499)',
    structure = SMILES('[CH2][CH]C([CH2])[CH2]'),
    E0 = (636.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1200.87],'cm^-1')),
        HinderedRotor(inertia=(0.112344,'amu*angstrom^2'), symmetry=1, barrier=(2.58301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112403,'amu*angstrom^2'), symmetry=1, barrier=(2.58438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252635,'amu*angstrom^2'), symmetry=1, barrier=(2.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252511,'amu*angstrom^2'), symmetry=1, barrier=(2.5843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95305,0.0434597,-3.02543e-05,1.30344e-08,-2.49244e-12,76589.1,26.1277], Tmin=(100,'K'), Tmax=(1189.32,'K')), NASAPolynomial(coeffs=[6.71646,0.027439,-1.00484e-05,1.70807e-09,-1.11583e-13,75456,2.32115], Tmin=(1189.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C([CH2])[CH2](646)',
    structure = SMILES('[CH2][CH][CH]C([CH2])[CH2]'),
    E0 = (806.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1304.76],'cm^-1')),
        HinderedRotor(inertia=(0.00319015,'amu*angstrom^2'), symmetry=1, barrier=(3.8543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318875,'amu*angstrom^2'), symmetry=1, barrier=(3.85345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0031889,'amu*angstrom^2'), symmetry=1, barrier=(3.85323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167613,'amu*angstrom^2'), symmetry=1, barrier=(3.85376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00319051,'amu*angstrom^2'), symmetry=1, barrier=(3.85453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52596,0.0581986,-6.27215e-05,4.88717e-08,-1.62492e-11,97126,32.6431], Tmin=(100,'K'), Tmax=(876.179,'K')), NASAPolynomial(coeffs=[4.1578,0.0383786,-1.54281e-05,2.72027e-09,-1.79882e-13,96964.3,22.0036], Tmin=(876.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(806.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][CH2](9040)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH2]'),
    E0 = (810.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1610.44,1610.45],'cm^-1')),
        HinderedRotor(inertia=(0.0851314,'amu*angstrom^2'), symmetry=1, barrier=(5.99533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0032575,'amu*angstrom^2'), symmetry=1, barrier=(5.99512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(6.49989e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0851266,'amu*angstrom^2'), symmetry=1, barrier=(5.99524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0336839,'amu*angstrom^2'), symmetry=1, barrier=(61.9916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75127,0.0531098,-4.22739e-05,2.36613e-08,-6.51376e-12,97552.1,31.6731], Tmin=(100,'K'), Tmax=(809.517,'K')), NASAPolynomial(coeffs=[4.44911,0.0397791,-1.75725e-05,3.3186e-09,-2.31342e-13,97115.3,19.2277], Tmin=(809.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(810.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][C]C([CH2])C([CH2])[CH2](9248)',
    structure = SMILES('[CH2][C]C([CH2])C([CH2])[CH2]'),
    E0 = (1045.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448502,0.0763222,-6.7345e-05,2.8234e-08,-2.3717e-12,125817,34.8779], Tmin=(100,'K'), Tmax=(830.826,'K')), NASAPolynomial(coeffs=[12.4302,0.0317596,-1.05829e-05,1.6986e-09,-1.07658e-13,123373,-23.4318], Tmin=(830.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1045.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C([CH2])[CH2](9249)',
    structure = SMILES('[CH]C([CH][CH2])C([CH2])[CH2]'),
    E0 = (1034.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200.087,1089.31,1317.82,1952.81,2130.09],'cm^-1')),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632966,0.0733946,-6.67761e-05,3.57233e-08,-8.02303e-12,124541,36.0823], Tmin=(100,'K'), Tmax=(1058.34,'K')), NASAPolynomial(coeffs=[10.4918,0.0361322,-1.39622e-05,2.45405e-09,-1.6401e-13,122455,-12.0396], Tmin=(1058.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1034.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C([CH2])[CH][CH2](9250)',
    structure = SMILES('[CH]C([CH2])C([CH2])[CH][CH2]'),
    E0 = (1034.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200.087,1089.31,1317.82,1952.81,2130.09],'cm^-1')),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632966,0.0733946,-6.67761e-05,3.57233e-08,-8.02303e-12,124541,36.7755], Tmin=(100,'K'), Tmax=(1058.34,'K')), NASAPolynomial(coeffs=[10.4918,0.0361322,-1.39622e-05,2.45405e-09,-1.6401e-13,122455,-11.3465], Tmin=(1058.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1034.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C([CH2])[CH2](9251)',
    structure = SMILES('[CH][CH]C([CH2])C([CH2])[CH2]'),
    E0 = (1034.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,198.428,996.551,1328.73,1829.52,2235.3],'cm^-1')),
        HinderedRotor(inertia=(0.13053,'amu*angstrom^2'), symmetry=1, barrier=(3.01505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13053,'amu*angstrom^2'), symmetry=1, barrier=(3.01505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13053,'amu*angstrom^2'), symmetry=1, barrier=(3.01505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13053,'amu*angstrom^2'), symmetry=1, barrier=(3.01505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13053,'amu*angstrom^2'), symmetry=1, barrier=(3.01505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13053,'amu*angstrom^2'), symmetry=1, barrier=(3.01505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.477867,0.075981,-7.46842e-05,4.37457e-08,-1.0583e-11,124528,36.6316], Tmin=(100,'K'), Tmax=(997.637,'K')), NASAPolynomial(coeffs=[10.8312,0.034469,-1.22678e-05,2.03562e-09,-1.30662e-13,122462,-13.2926], Tmin=(997.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1034.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])C1CC1[CH2](5208)',
    structure = SMILES('[CH2]C([CH2])C1CC1[CH2]'),
    E0 = (535.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,669.3,784.116,784.116,784.116,784.116,784.116,784.116,784.116,784.116,784.116,784.116,2517.52],'cm^-1')),
        HinderedRotor(inertia=(0.00119983,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119983,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119983,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119983,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04271,0.0496222,2.37629e-05,-7.11488e-08,3.35016e-11,64475.4,29.6271], Tmin=(100,'K'), Tmax=(907.029,'K')), NASAPolynomial(coeffs=[15.3669,0.0272927,-6.84888e-06,9.922e-10,-6.46684e-14,60197,-47.3417], Tmin=(907.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])C1CC1(9252)',
    structure = SMILES('[CH2][CH]C([CH2])C1CC1'),
    E0 = (537.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14361,0.0492651,1.31171e-05,-4.78936e-08,2.08621e-11,64711.5,31.7081], Tmin=(100,'K'), Tmax=(986.826,'K')), NASAPolynomial(coeffs=[12.7377,0.0342208,-1.25814e-05,2.27721e-09,-1.59931e-13,60867.5,-31.9553], Tmin=(986.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1CC([CH2])C1[CH2](5209)',
    structure = SMILES('[CH2]C1CC([CH2])C1[CH2]'),
    E0 = (530.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32323,0.0406985,4.96372e-05,-9.63123e-08,4.17599e-11,63936.3,28.6336], Tmin=(100,'K'), Tmax=(916.853,'K')), NASAPolynomial(coeffs=[14.9238,0.0284402,-7.32797e-06,1.11184e-09,-7.53851e-14,59463.6,-46.5915], Tmin=(916.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CCC1[CH2](5275)',
    structure = SMILES('[CH2][CH]C1CCC1[CH2]'),
    E0 = (529.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44015,0.0401354,3.98002e-05,-7.42515e-08,2.97041e-11,63769,29.9649], Tmin=(100,'K'), Tmax=(982.689,'K')), NASAPolynomial(coeffs=[12.2781,0.0353991,-1.30795e-05,2.40153e-09,-1.71048e-13,59737.6,-31.8069], Tmin=(982.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(=C)C([CH2])[CH2](9253)',
    structure = SMILES('[CH2]CC(=C)C([CH2])[CH2]'),
    E0 = (509.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539552,0.0689308,-4.23226e-05,6.63737e-09,2.89082e-12,61424.2,32.823], Tmin=(100,'K'), Tmax=(968.886,'K')), NASAPolynomial(coeffs=[12.7498,0.0338736,-1.18164e-05,2.00118e-09,-1.32792e-13,58337.5,-29.4174], Tmin=(968.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(C)C([CH2])[CH2](4089)',
    structure = SMILES('[CH2]C([CH2])[C](C)C=C'),
    E0 = (440.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,306.003,1441.35],'cm^-1')),
        HinderedRotor(inertia=(8.11663e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126145,'amu*angstrom^2'), symmetry=1, barrier=(8.40394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00569836,'amu*angstrom^2'), symmetry=1, barrier=(8.40097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03113,'amu*angstrom^2'), symmetry=1, barrier=(68.9845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03761,'amu*angstrom^2'), symmetry=1, barrier=(68.99,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904511,0.0590154,-1.4749e-05,-2.12076e-08,1.27911e-11,53113.6,28.5207], Tmin=(100,'K'), Tmax=(938.632,'K')), NASAPolynomial(coeffs=[11.8409,0.0347462,-1.16606e-05,1.95166e-09,-1.29852e-13,50076.6,-28.79], Tmin=(938.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]CC([CH2])C([CH2])=C(9254)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])=C'),
    E0 = (454.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439407,0.0686603,-3.61972e-05,-9.86269e-10,5.35541e-12,54848.8,31.0138], Tmin=(100,'K'), Tmax=(1007.43,'K')), NASAPolynomial(coeffs=[14.2461,0.0331702,-1.21347e-05,2.13587e-09,-1.45657e-13,51086.1,-40.5656], Tmin=(1007.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])=C(4090)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])=C'),
    E0 = (444.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0125055,'amu*angstrom^2'), symmetry=1, barrier=(22.7494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880852,'amu*angstrom^2'), symmetry=1, barrier=(2.02525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868833,'amu*angstrom^2'), symmetry=1, barrier=(1.99762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125217,'amu*angstrom^2'), symmetry=1, barrier=(22.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125178,'amu*angstrom^2'), symmetry=1, barrier=(22.7631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.523855,0.066626,-3.72749e-05,6.24285e-09,9.963e-13,53577.7,31.4961], Tmin=(100,'K'), Tmax=(1193.41,'K')), NASAPolynomial(coeffs=[13.9008,0.0349916,-1.41067e-05,2.56982e-09,-1.76018e-13,49444.8,-39.344], Tmin=(1193.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cs_S) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])[CH2](5212)',
    structure = SMILES('[CH2]C=C([CH2])C([CH2])[CH2]'),
    E0 = (593.992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,677.293],'cm^-1')),
        HinderedRotor(inertia=(0.263955,'amu*angstrom^2'), symmetry=1, barrier=(84.0881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0085887,'amu*angstrom^2'), symmetry=1, barrier=(84.0126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0786157,'amu*angstrom^2'), symmetry=1, barrier=(1.80753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258323,'amu*angstrom^2'), symmetry=1, barrier=(83.8601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.64165,'amu*angstrom^2'), symmetry=1, barrier=(84.1452,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.537746,0.0664897,-3.12729e-05,-1.07077e-08,1.07207e-11,71574.3,30.499], Tmin=(100,'K'), Tmax=(931.076,'K')), NASAPolynomial(coeffs=[15.0887,0.0283383,-9.05616e-06,1.48607e-09,-9.87628e-14,67808.8,-44.3323], Tmin=(931.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]([CH2])C([CH2])C=C(5213)',
    structure = SMILES('[CH2][C]([CH2])C([CH2])C=C'),
    E0 = (699.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,180,1674.9],'cm^-1')),
        HinderedRotor(inertia=(0.183777,'amu*angstrom^2'), symmetry=1, barrier=(4.22539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00212315,'amu*angstrom^2'), symmetry=1, barrier=(4.22622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183929,'amu*angstrom^2'), symmetry=1, barrier=(4.2289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0320922,'amu*angstrom^2'), symmetry=1, barrier=(63.8773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183739,'amu*angstrom^2'), symmetry=1, barrier=(4.22451,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759497,0.0731758,-7.83801e-05,5.44667e-08,-1.58121e-11,84202,33.7049], Tmin=(100,'K'), Tmax=(926.11,'K')), NASAPolynomial(coeffs=[7.55971,0.0387103,-1.43057e-05,2.4026e-09,-1.54131e-13,83161,2.59928], Tmin=(926.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(699.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=CC([CH2])[CH2](644)',
    structure = SMILES('[CH2]C=CC([CH2])[CH2]'),
    E0 = (477.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,397.773],'cm^-1')),
        HinderedRotor(inertia=(0.00475792,'amu*angstrom^2'), symmetry=1, barrier=(0.534766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607471,'amu*angstrom^2'), symmetry=1, barrier=(72.6648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0695163,'amu*angstrom^2'), symmetry=1, barrier=(7.65916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0302878,'amu*angstrom^2'), symmetry=1, barrier=(72.6694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3691,0.049131,-9.81136e-06,-2.36824e-08,1.37799e-11,57516.4,26.478], Tmin=(100,'K'), Tmax=(924.208,'K')), NASAPolynomial(coeffs=[11.9922,0.0255947,-8.03361e-06,1.30776e-09,-8.6727e-14,54594.4,-29.1203], Tmin=(924.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C=C(496)',
    structure = SMILES('[CH2]C([CH2])C=C'),
    E0 = (361.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,444.303],'cm^-1')),
        HinderedRotor(inertia=(0.0434666,'amu*angstrom^2'), symmetry=1, barrier=(6.11473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.043652,'amu*angstrom^2'), symmetry=1, barrier=(6.1174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0394383,'amu*angstrom^2'), symmetry=1, barrier=(69.5439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03654,0.0357961,3.13511e-06,-3.00403e-08,1.51068e-11,43603,21.9965], Tmin=(100,'K'), Tmax=(906.7,'K')), NASAPolynomial(coeffs=[9.64737,0.0219243,-6.51387e-06,1.02236e-09,-6.65462e-14,41412.9,-18.4424], Tmin=(906.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2][CH][CH2](497)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.000863049,'amu*angstrom^2'), symmetry=1, barrier=(2.40754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000864751,'amu*angstrom^2'), symmetry=1, barrier=(2.41365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42318,0.0135015,1.25385e-06,-4.52647e-09,1.21417e-12,58330.2,15.4092], Tmin=(100,'K'), Tmax=(1578.37,'K')), NASAPolynomial(coeffs=[5.16145,0.0152713,-6.29635e-06,1.14122e-09,-7.61383e-14,57012.3,3.79304], Tmin=(1578.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH2](5531)',
    structure = SMILES('[CH2][CH][CH][CH2]'),
    E0 = (655.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1802.64],'cm^-1')),
        HinderedRotor(inertia=(0.00215831,'amu*angstrom^2'), symmetry=1, barrier=(4.96293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00215621,'amu*angstrom^2'), symmetry=1, barrier=(4.95965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214837,'amu*angstrom^2'), symmetry=1, barrier=(4.95028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82147,0.0344596,-5.62911e-05,6.39381e-08,-2.62695e-11,78855.1,21.3291], Tmin=(100,'K'), Tmax=(865.068,'K')), NASAPolynomial(coeffs=[-1.36886,0.0321046,-1.45274e-05,2.71435e-09,-1.84229e-13,80393.2,45.6373], Tmin=(865.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C([CH2])[CH2](9255)',
    structure = SMILES('[CH2][CH][C]([CH2])C([CH2])[CH2]'),
    E0 = (976.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,1132.57,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716368,0.0803268,-0.000109836,9.48947e-08,-3.22458e-11,117589,37.6941], Tmin=(100,'K'), Tmax=(904.093,'K')), NASAPolynomial(coeffs=[4.28057,0.0448227,-1.81877e-05,3.17085e-09,-2.06047e-13,117751,25.3187], Tmin=(904.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(976.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]([CH2])[CH2](9256)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])[CH2]'),
    E0 = (976.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,1132.57,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0103496,'amu*angstrom^2'), symmetry=1, barrier=(5.94567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716368,0.0803268,-0.000109836,9.48947e-08,-3.22458e-11,117589,37.6941], Tmin=(100,'K'), Tmax=(904.093,'K')), NASAPolynomial(coeffs=[4.28057,0.0448227,-1.81877e-05,3.17085e-09,-2.06047e-13,117751,25.3187], Tmin=(904.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(976.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])C([CH2])[CH2](5204)',
    structure = SMILES('[CH2]C[C]([CH2])C([CH2])[CH2]'),
    E0 = (782.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,878.673,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510331,0.081721,-9.90187e-05,7.73506e-08,-2.45549e-11,94201.4,35.8774], Tmin=(100,'K'), Tmax=(912.514,'K')), NASAPolynomial(coeffs=[6.39508,0.0438586,-1.69449e-05,2.89768e-09,-1.87012e-13,93629.8,10.7784], Tmin=(912.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C](C)C([CH2])[CH2](4091)',
    structure = SMILES('[CH2][CH][C](C)C([CH2])[CH2]'),
    E0 = (771.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,949.828,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754627,0.0779817,-9.51037e-05,7.97588e-08,-2.76948e-11,92923.3,35.7755], Tmin=(100,'K'), Tmax=(858.426,'K')), NASAPolynomial(coeffs=[4.13625,0.0488031,-2.06656e-05,3.73589e-09,-2.50362e-13,92837.3,22.8576], Tmin=(858.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]([CH2])C(9065)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])C'),
    E0 = (771.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,949.828,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754627,0.0779817,-9.51037e-05,7.97588e-08,-2.76948e-11,92923.3,36.4687], Tmin=(100,'K'), Tmax=(858.426,'K')), NASAPolynomial(coeffs=[4.13625,0.0488031,-2.06656e-05,3.73589e-09,-2.50362e-13,92837.3,23.5508], Tmin=(858.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[C]([CH2])[CH2](5205)',
    structure = SMILES('[CH2]CC([CH2])[C]([CH2])[CH2]'),
    E0 = (782.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,878.673,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156122,'amu*angstrom^2'), symmetry=1, barrier=(6.78638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510331,0.081721,-9.90187e-05,7.73506e-08,-2.45549e-11,94201.4,35.8774], Tmin=(100,'K'), Tmax=(912.514,'K')), NASAPolynomial(coeffs=[6.39508,0.0438586,-1.69449e-05,2.89768e-09,-1.87012e-13,93629.8,10.7784], Tmin=(912.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)[C]([CH2])[CH2](4092)',
    structure = SMILES('[CH2][CH]C(C)[C]([CH2])[CH2]'),
    E0 = (771.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,949.828,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754627,0.0779817,-9.51037e-05,7.97588e-08,-2.76948e-11,92923.3,35.7755], Tmin=(100,'K'), Tmax=(858.426,'K')), NASAPolynomial(coeffs=[4.13625,0.0488031,-2.06656e-05,3.73589e-09,-2.50362e-13,92837.3,22.8576], Tmin=(858.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C([CH2])C(9066)',
    structure = SMILES('[CH2][CH][C]([CH2])C([CH2])C'),
    E0 = (771.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,949.828,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131786,'amu*angstrom^2'), symmetry=1, barrier=(6.30775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754627,0.0779817,-9.51037e-05,7.97588e-08,-2.76948e-11,92923.3,36.4687], Tmin=(100,'K'), Tmax=(858.426,'K')), NASAPolynomial(coeffs=[4.13625,0.0488031,-2.06656e-05,3.73589e-09,-2.50362e-13,92837.3,23.5508], Tmin=(858.426,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C([CH2])[CH2](4123)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])[CH2]'),
    E0 = (771.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,941.124,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686567,0.0795013,-9.90437e-05,8.21785e-08,-2.76038e-11,92906.1,36.0148], Tmin=(100,'K'), Tmax=(895.828,'K')), NASAPolynomial(coeffs=[4.46691,0.0471529,-1.89776e-05,3.3187e-09,-2.171e-13,92849.5,21.657], Tmin=(895.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C([CH2])[CH]C(4124)',
    structure = SMILES('[CH2][C]([CH2])C([CH2])[CH]C'),
    E0 = (771.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,941.124,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131647,'amu*angstrom^2'), symmetry=1, barrier=(6.29142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686567,0.0795013,-9.90437e-05,8.21785e-08,-2.76038e-11,92906.1,36.0148], Tmin=(100,'K'), Tmax=(895.828,'K')), NASAPolynomial(coeffs=[4.46691,0.0471529,-1.89776e-05,3.3187e-09,-2.171e-13,92849.5,21.657], Tmin=(895.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (791.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (791.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (948.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (948.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1249.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1244.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1248.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1256.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1246.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1246.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1246.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (796.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (799.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (799.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (799.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (814.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (814.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (814.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (854.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (854.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (854.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (813.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (922.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (898.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (919.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (915.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (796.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (832.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1140.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1188.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (1188.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (949.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (932.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (932.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (907.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (909.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (909.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (948.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (922.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['allyl(82)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]C([CH2])C([CH2])C=C(4122)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]CC([CH2])[CH][CH2](4108)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C([CH2])[CH2](499)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C([CH2])[CH2](646)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])[CH][CH2](9040)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C]C([CH2])C([CH2])[CH2](9248)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([CH][CH2])C([CH2])[CH2](9249)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([CH2])C([CH2])[CH][CH2](9250)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH]C([CH2])C([CH2])[CH2](9251)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]C([CH2])C1CC1[CH2](5208)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]C([CH2])C1CC1(9252)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]C1CC([CH2])C1[CH2](5209)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]C1CCC1[CH2](5275)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]CC(=C)C([CH2])[CH2](9253)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]C=C(C)C([CH2])[CH2](4089)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][C](C)C([CH2])C=C(4991)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][C](C=C)C([CH2])C(4992)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]CC([CH2])C([CH2])=C(9254)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]C(C)C([CH2])=C(4090)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C](C=C)C([CH2])[CH2](5212)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(20)', '[CH2]C=CC([CH2])[CH2](644)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH][CH2](502)', '[CH2]C([CH2])C=C(496)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])C=C(4950)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['allyl(82)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.0036002,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
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
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C([CH2])[CH2](9255)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[C]([CH2])[CH2](9256)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]C[C]([CH2])C([CH2])[CH2](5204)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH][C](C)C([CH2])[CH2](4091)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]C([CH2])[C]([CH2])C(9065)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2]CC([CH2])[C]([CH2])[CH2](5205)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]C(C)[C]([CH2])[CH2](4092)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH][C]([CH2])C([CH2])C(9066)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][C]([CH]C)C([CH2])[CH2](4123)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][C]([CH2])C([CH2])[CH]C(4124)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2053',
    isomers = [
        '[CH2][CH]C([CH2])C([CH2])[CH2](4093)',
    ],
    reactants = [
        ('allyl(82)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2053',
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


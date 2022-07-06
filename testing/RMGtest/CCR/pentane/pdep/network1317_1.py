species(
    label = '[CH2][CH]CC([CH2])C=C(4148)',
    structure = SMILES('[CH2][CH]CC([CH2])C=C'),
    E0 = (508.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,268.431,930.189,2053.03],'cm^-1')),
        HinderedRotor(inertia=(0.0810729,'amu*angstrom^2'), symmetry=1, barrier=(3.35781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810729,'amu*angstrom^2'), symmetry=1, barrier=(3.35781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810729,'amu*angstrom^2'), symmetry=1, barrier=(3.35781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810729,'amu*angstrom^2'), symmetry=1, barrier=(3.35781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810729,'amu*angstrom^2'), symmetry=1, barrier=(3.35781,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.960089,0.0647244,-4.26874e-05,1.5814e-08,-2.54031e-12,61322.3,33.7778], Tmin=(100,'K'), Tmax=(1387.17,'K')), NASAPolynomial(coeffs=[9.69163,0.0395464,-1.54614e-05,2.72936e-09,-1.82155e-13,58899.9,-11.2044], Tmin=(1387.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C(C=C)CC=C(4875)',
    structure = SMILES('[CH2]C(C=C)CC=C'),
    E0 = (237.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,625.486,651.228,657.699],'cm^-1')),
        HinderedRotor(inertia=(0.437694,'amu*angstrom^2'), symmetry=1, barrier=(10.0634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168645,'amu*angstrom^2'), symmetry=1, barrier=(3.87748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433536,'amu*angstrom^2'), symmetry=1, barrier=(9.96785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807781,'amu*angstrom^2'), symmetry=1, barrier=(24.5586,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.834958,0.0582886,-1.11822e-05,-2.41059e-08,1.28659e-11,28632.1,30.1822], Tmin=(100,'K'), Tmax=(994.87,'K')), NASAPolynomial(coeffs=[13.1761,0.0338944,-1.24346e-05,2.21896e-09,-1.53504e-13,24928.2,-35.5671], Tmin=(994.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC[CH]C=C(4160)',
    structure = SMILES('[CH2][CH]CCC=C[CH2]'),
    E0 = (451.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948154,0.0638066,-3.81353e-05,1.15283e-08,-1.44607e-12,54364.8,32.1312], Tmin=(100,'K'), Tmax=(1753.37,'K')), NASAPolynomial(coeffs=[13.0017,0.0363086,-1.46109e-05,2.58389e-09,-1.70758e-13,50138,-32.789], Tmin=(1753.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C[CH]CC=C(5244)',
    structure = SMILES('[CH2][CH]C[CH]CC=C'),
    E0 = (506.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.419,0.0530823,-2.41139e-05,4.6519e-09,-3.32175e-13,60929.5,29.458], Tmin=(100,'K'), Tmax=(2959.84,'K')), NASAPolynomial(coeffs=[39.6324,0.00734902,-3.2468e-06,4.72089e-10,-2.30743e-14,36904,-193.829], Tmin=(2959.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][C]CC([CH2])C=C(5245)',
    structure = SMILES('[CH2][C]CC([CH2])C=C'),
    E0 = (762.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355835,0.0716049,-5.6303e-05,2.37928e-08,-4.09487e-12,91871,33.1256], Tmin=(100,'K'), Tmax=(1380.1,'K')), NASAPolynomial(coeffs=[14.7669,0.0298363,-1.09051e-05,1.86283e-09,-1.22284e-13,87893.3,-41.042], Tmin=(1380.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C)C[CH][CH2](5246)',
    structure = SMILES('[CH]C(C=C)C[CH][CH2]'),
    E0 = (752.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00253,0.0647706,-4.508e-05,1.69811e-08,-2.72666e-12,90562,32.8473], Tmin=(100,'K'), Tmax=(1398.97,'K')), NASAPolynomial(coeffs=[10.7225,0.0369784,-1.52804e-05,2.78021e-09,-1.88898e-13,87842.5,-17.3093], Tmin=(1398.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC([CH2])C=C(5247)',
    structure = SMILES('[CH][CH]CC([CH2])C=C'),
    E0 = (751.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884291,0.0668622,-5.10869e-05,2.24695e-08,-4.23142e-12,90547.3,33.2699], Tmin=(100,'K'), Tmax=(1225.24,'K')), NASAPolynomial(coeffs=[10.0292,0.0370073,-1.45369e-05,2.58226e-09,-1.73588e-13,88306.3,-12.7065], Tmin=(1225.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC(C=C)C1(5237)',
    structure = SMILES('[CH2]C1CC(C=C)C1'),
    E0 = (258.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52854,0.0324583,7.30303e-05,-1.16822e-07,4.69522e-11,31185.2,25.8132], Tmin=(100,'K'), Tmax=(951.958,'K')), NASAPolynomial(coeffs=[15.2651,0.0297816,-9.48233e-06,1.70042e-09,-1.24695e-13,26075.9,-52.8808], Tmin=(951.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCC(=C)C=C(5248)',
    structure = SMILES('[CH2]CCC(=C)C=C'),
    E0 = (215.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464399,0.0656976,-2.47619e-05,-1.44776e-08,1.02152e-11,26014.6,28.2297], Tmin=(100,'K'), Tmax=(1008.97,'K')), NASAPolynomial(coeffs=[15.6299,0.0312515,-1.17244e-05,2.12971e-09,-1.49091e-13,21647.3,-51.5471], Tmin=(1008.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC(C)C=C(5249)',
    structure = SMILES('[CH2]C=CC(C)C=C'),
    E0 = (170.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51279,0.067028,-3.84304e-05,7.5765e-09,5.44575e-13,20658.4,27.0237], Tmin=(100,'K'), Tmax=(1202.28,'K')), NASAPolynomial(coeffs=[13.7648,0.035208,-1.4039e-05,2.53977e-09,-1.73151e-13,16585.1,-43.0391], Tmin=(1202.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])C[CH2](5250)',
    structure = SMILES('[CH2][CH]C[C]([CH2])C[CH2]'),
    E0 = (777.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,230.801,3157.68,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0125513,'amu*angstrom^2'), symmetry=1, barrier=(2.25359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125513,'amu*angstrom^2'), symmetry=1, barrier=(2.25359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125513,'amu*angstrom^2'), symmetry=1, barrier=(2.25359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125513,'amu*angstrom^2'), symmetry=1, barrier=(2.25359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125513,'amu*angstrom^2'), symmetry=1, barrier=(2.25359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0125513,'amu*angstrom^2'), symmetry=1, barrier=(2.25359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.774016,0.0809608,-0.000111392,1.03073e-07,-3.78913e-11,93612.4,36.8614], Tmin=(100,'K'), Tmax=(861.187,'K')), NASAPolynomial(coeffs=[1.79688,0.0533831,-2.35985e-05,4.33118e-09,-2.91687e-13,94282.7,36.9941], Tmin=(861.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(777.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(RCCJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])C[CH2](5251)',
    structure = SMILES('[CH2][CH][CH]C([CH2])C[CH2]'),
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
    label = '[CH2][CH][CH]C([CH2])[CH]C(4150)',
    structure = SMILES('[CH2][CH][CH]C([CH2])[CH]C'),
    E0 = (775.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1455.87,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10013,0.069552,-7.39774e-05,6.16456e-08,-2.28945e-11,93414.1,37.884], Tmin=(100,'K'), Tmax=(797.472,'K')), NASAPolynomial(coeffs=[2.76063,0.0518444,-2.30295e-05,4.3071e-09,-2.9626e-13,93447.5,32.1186], Tmin=(797.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(775.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC1[CH]CC1(5252)',
    structure = SMILES('[CH2][CH]CC1[CH]CC1'),
    E0 = (521.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5491,0.0432809,1.53255e-05,-3.75513e-08,1.38934e-11,62761.6,31.2982], Tmin=(100,'K'), Tmax=(1089.31,'K')), NASAPolynomial(coeffs=[8.4488,0.0425307,-1.74966e-05,3.25576e-09,-2.27183e-13,59799.7,-9.27453], Tmin=(1089.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])C1(5227)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])C1'),
    E0 = (447.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79424,0.0295658,7.13455e-05,-1.10585e-07,4.43356e-11,53905.3,28.8813], Tmin=(100,'K'), Tmax=(938.097,'K')), NASAPolynomial(coeffs=[12.4009,0.0323445,-9.85695e-06,1.67205e-09,-1.17752e-13,49803,-32.8704], Tmin=(938.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]CC[CH]C1(5253)',
    structure = SMILES('[CH2]C1[CH]CC[CH]C1'),
    E0 = (407.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9442,0.0325415,4.73336e-05,-6.57075e-08,2.20629e-11,49041.4,25.1849], Tmin=(100,'K'), Tmax=(1093.08,'K')), NASAPolynomial(coeffs=[6.58039,0.0490039,-2.11294e-05,4.0251e-09,-2.84479e-13,46030.8,-6.72957], Tmin=(1093.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclohexane) + radical(Cs_S) + radical(cyclohexane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC1CC1[CH2](5254)',
    structure = SMILES('[CH2][CH]CC1CC1[CH2]'),
    E0 = (533.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18626,0.0502917,6.39475e-06,-3.91097e-08,1.76173e-11,64294.1,31.0406], Tmin=(100,'K'), Tmax=(979.839,'K')), NASAPolynomial(coeffs=[11.4048,0.0357259,-1.28691e-05,2.27536e-09,-1.5683e-13,60988.3,-24.7006], Tmin=(979.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C1C[CH]CC1[CH2](5255)',
    structure = SMILES('[CH2]C1C[CH]CC1[CH2]'),
    E0 = (438.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87887,0.0309914,5.98151e-05,-9.4117e-08,3.76557e-11,52857.5,27.8925], Tmin=(100,'K'), Tmax=(935.171,'K')), NASAPolynomial(coeffs=[10.0526,0.0356735,-1.12825e-05,1.89767e-09,-1.30219e-13,49595.3,-20.2615], Tmin=(935.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(cyclopentane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC(=C)C=C(5256)',
    structure = SMILES('[CH2][CH]CC(=C)C=C'),
    E0 = (409.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3000,3100,440,815,1455,1000,208.162,208.162,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00389045,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610721,'amu*angstrom^2'), symmetry=1, barrier=(18.7789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610721,'amu*angstrom^2'), symmetry=1, barrier=(18.7789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.042029,'amu*angstrom^2'), symmetry=1, barrier=(47.9429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519406,0.0674144,-4.86029e-05,1.85157e-08,-2.88169e-12,49396.1,30.7796], Tmin=(100,'K'), Tmax=(1506.62,'K')), NASAPolynomial(coeffs=[14.7235,0.0297028,-1.10566e-05,1.90162e-09,-1.24805e-13,45116.1,-43.5686], Tmin=(1506.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])C=C(5257)',
    structure = SMILES('[CH2]C=CC([CH2])C=C'),
    E0 = (375.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,3850.91],'cm^-1')),
        HinderedRotor(inertia=(0.530733,'amu*angstrom^2'), symmetry=1, barrier=(12.2026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07697,'amu*angstrom^2'), symmetry=1, barrier=(24.7617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743207,'amu*angstrom^2'), symmetry=1, barrier=(89.5211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.89708,'amu*angstrom^2'), symmetry=1, barrier=(89.6014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547792,0.0685847,-5.08023e-05,2.03371e-08,-3.35663e-12,45320.9,28.6747], Tmin=(100,'K'), Tmax=(1419.13,'K')), NASAPolynomial(coeffs=[13.6133,0.0317576,-1.18766e-05,2.05085e-09,-1.3524e-13,41612.6,-38.9326], Tmin=(1419.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = '[CH2][CH]CC=C(508)',
    structure = SMILES('[CH2][CH]CC=C'),
    E0 = (359.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1243.57],'cm^-1')),
        HinderedRotor(inertia=(0.00331557,'amu*angstrom^2'), symmetry=1, barrier=(3.65555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160798,'amu*angstrom^2'), symmetry=1, barrier=(3.69707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161411,'amu*angstrom^2'), symmetry=1, barrier=(3.71116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13435,0.0369817,-1.51234e-05,1.33725e-09,3.77008e-13,43295.2,23.5599], Tmin=(100,'K'), Tmax=(1532.36,'K')), NASAPolynomial(coeffs=[8.68459,0.025936,-1.02358e-05,1.78847e-09,-1.17139e-13,40577.1,-13.1557], Tmin=(1532.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH][CH2](510)',
    structure = SMILES('[CH2][CH]C[CH][CH2]'),
    E0 = (631.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1407.99,1417.06],'cm^-1')),
        HinderedRotor(inertia=(0.00560336,'amu*angstrom^2'), symmetry=1, barrier=(7.95035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00563889,'amu*angstrom^2'), symmetry=1, barrier=(7.95031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00557143,'amu*angstrom^2'), symmetry=1, barrier=(7.9482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00562946,'amu*angstrom^2'), symmetry=1, barrier=(7.93664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26906,0.0435841,-4.87505e-05,4.66768e-08,-1.86942e-11,75984.9,26.4118], Tmin=(100,'K'), Tmax=(836.615,'K')), NASAPolynomial(coeffs=[0.578549,0.0388266,-1.7199e-05,3.18971e-09,-2.17339e-13,76717.2,36.9515], Tmin=(836.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])C=C(5258)',
    structure = SMILES('[CH2][CH]CC([CH2])=C[CH2]'),
    E0 = (587.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1499.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.195651,'amu*angstrom^2'), symmetry=1, barrier=(8.76227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704221,'amu*angstrom^2'), symmetry=1, barrier=(31.5379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704202,'amu*angstrom^2'), symmetry=1, barrier=(31.5379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0715767,'amu*angstrom^2'), symmetry=1, barrier=(114.224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0715766,'amu*angstrom^2'), symmetry=1, barrier=(114.224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756853,0.0664382,-4.5934e-05,1.67474e-08,-2.54131e-12,70757.4,31.612], Tmin=(100,'K'), Tmax=(1505.17,'K')), NASAPolynomial(coeffs=[12.967,0.0339898,-1.35972e-05,2.42502e-09,-1.62456e-13,67081.7,-32.2879], Tmin=(1505.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])C=C(5259)',
    structure = SMILES('[CH2][CH][CH]C([CH2])C=C'),
    E0 = (703.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,275.369,1343.59,2995.74],'cm^-1')),
        HinderedRotor(inertia=(0.00249673,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00249673,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00249673,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00249673,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00249673,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39061,0.0605369,-4.31916e-05,1.94384e-08,-4.12707e-12,84700.7,34.8008], Tmin=(100,'K'), Tmax=(1029.57,'K')), NASAPolynomial(coeffs=[5.764,0.0435456,-1.84363e-05,3.40874e-09,-2.34707e-13,83800.2,13.5743], Tmin=(1029.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(703.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])[C]=C(5260)',
    structure = SMILES('[CH2][CH]CC([CH2])[C]=C'),
    E0 = (746.782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,284.834,962.061,2120.18],'cm^-1')),
        HinderedRotor(inertia=(0.0708138,'amu*angstrom^2'), symmetry=1, barrier=(3.30232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0708138,'amu*angstrom^2'), symmetry=1, barrier=(3.30232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0708138,'amu*angstrom^2'), symmetry=1, barrier=(3.30232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0708138,'amu*angstrom^2'), symmetry=1, barrier=(3.30232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0708138,'amu*angstrom^2'), symmetry=1, barrier=(3.30232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05735,0.0676219,-6.09214e-05,3.60595e-08,-9.5482e-12,89920.6,33.9159], Tmin=(100,'K'), Tmax=(880.574,'K')), NASAPolynomial(coeffs=[6.66682,0.0421405,-1.75148e-05,3.19664e-09,-2.18113e-13,88932.7,7.56696], Tmin=(880.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(746.782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C[CH][CH2](5261)',
    structure = SMILES('[CH]=CC([CH2])C[CH][CH2]'),
    E0 = (756.037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,180,2319.45],'cm^-1')),
        HinderedRotor(inertia=(0.00317465,'amu*angstrom^2'), symmetry=1, barrier=(12.1234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.527235,'amu*angstrom^2'), symmetry=1, barrier=(12.1222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0196179,'amu*angstrom^2'), symmetry=1, barrier=(74.9083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147961,'amu*angstrom^2'), symmetry=1, barrier=(74.9112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124505,'amu*angstrom^2'), symmetry=1, barrier=(2.86261,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995263,0.066887,-5.43688e-05,2.6653e-08,-5.69982e-12,91037.7,34.0617], Tmin=(100,'K'), Tmax=(1080.97,'K')), NASAPolynomial(coeffs=[8.51479,0.0390623,-1.57583e-05,2.8411e-09,-1.92807e-13,89412,-2.80123], Tmin=(1080.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(756.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]CC(C)=C[CH2](4104)',
    structure = SMILES('[CH2][CH]C[C](C)C=C'),
    E0 = (435.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,281.543,955.88,2102.45],'cm^-1')),
        HinderedRotor(inertia=(0.07278,'amu*angstrom^2'), symmetry=1, barrier=(3.31603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.07278,'amu*angstrom^2'), symmetry=1, barrier=(3.31603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.07278,'amu*angstrom^2'), symmetry=1, barrier=(3.31603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.07278,'amu*angstrom^2'), symmetry=1, barrier=(3.31603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.07278,'amu*angstrom^2'), symmetry=1, barrier=(3.31603,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948327,0.0604438,-3.27883e-05,8.80954e-09,-9.63893e-13,52534.7,30.3197], Tmin=(100,'K'), Tmax=(2024.99,'K')), NASAPolynomial(coeffs=[14.9777,0.032732,-1.22614e-05,2.05185e-09,-1.29627e-13,46852.7,-47.263], Tmin=(2024.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C([CH2])C=C(5262)',
    structure = SMILES('[CH2]C[CH]C([CH2])C=C'),
    E0 = (509.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,244.808,884.016,1958.47],'cm^-1')),
        HinderedRotor(inertia=(0.097905,'amu*angstrom^2'), symmetry=1, barrier=(3.37267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097905,'amu*angstrom^2'), symmetry=1, barrier=(3.37267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097905,'amu*angstrom^2'), symmetry=1, barrier=(3.37267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097905,'amu*angstrom^2'), symmetry=1, barrier=(3.37267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097905,'amu*angstrom^2'), symmetry=1, barrier=(3.37267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729455,0.0657575,-4.24415e-05,1.44274e-08,-2.04309e-12,61345.5,34.437], Tmin=(100,'K'), Tmax=(1599,'K')), NASAPolynomial(coeffs=[13.0727,0.0348803,-1.34761e-05,2.35095e-09,-1.54971e-13,57398.2,-30.9055], Tmin=(1599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C=C[CH2](4105)',
    structure = SMILES('[CH2][CH]C(C)C=C[CH2]'),
    E0 = (443.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1226.02],'cm^-1')),
        HinderedRotor(inertia=(0.182504,'amu*angstrom^2'), symmetry=1, barrier=(4.19614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218248,'amu*angstrom^2'), symmetry=1, barrier=(23.2796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00393486,'amu*angstrom^2'), symmetry=1, barrier=(4.19745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182714,'amu*angstrom^2'), symmetry=1, barrier=(4.20096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01233,'amu*angstrom^2'), symmetry=1, barrier=(23.2754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722107,0.0622069,-2.74146e-05,-2.0837e-09,3.45469e-12,53431.8,32.5258], Tmin=(100,'K'), Tmax=(1157.75,'K')), NASAPolynomial(coeffs=[12.8894,0.0363317,-1.4831e-05,2.72869e-09,-1.88342e-13,49531.2,-32.6348], Tmin=(1157.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC(C)[C]=C(5263)',
    structure = SMILES('[CH2][CH]CC(C)[C]=C'),
    E0 = (541.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25766,0.063577,-4.09812e-05,1.4711e-08,-2.37399e-12,65247.1,31.3991], Tmin=(100,'K'), Tmax=(1314.44,'K')), NASAPolynomial(coeffs=[7.62959,0.0441865,-1.88533e-05,3.48809e-09,-2.39448e-13,63572,-1.08399], Tmin=(1314.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC[C]([CH2])C=C(5264)',
    structure = SMILES('[CH2]C=C([CH2])CC[CH2]'),
    E0 = (392.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.405852,0.0680415,-3.29821e-05,-3.03876e-09,5.27786e-12,47389.2,30.1364], Tmin=(100,'K'), Tmax=(1066.22,'K')), NASAPolynomial(coeffs=[14.8757,0.0339317,-1.3378e-05,2.45035e-09,-1.70341e-13,43156.8,-45.9779], Tmin=(1066.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH][CH]C)C=C(5043)',
    structure = SMILES('[CH2]C([CH][CH]C)C=C'),
    E0 = (498.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,402.67,1554.62,3088.64],'cm^-1')),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34038,0.0602115,-3.50779e-05,1.1048e-08,-1.54084e-12,60017.7,33.1737], Tmin=(100,'K'), Tmax=(1488.46,'K')), NASAPolynomial(coeffs=[7.87184,0.0426592,-1.73893e-05,3.12547e-09,-2.10162e-13,58073.3,-0.93461], Tmin=(1488.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCC([CH2])[C]=C(5265)',
    structure = SMILES('[CH2]CCC([CH2])[C]=C'),
    E0 = (552.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.630078,0.070238,-5.17878e-05,2.11694e-08,-3.64498e-12,66555.1,32.7032], Tmin=(100,'K'), Tmax=(1340.94,'K')), NASAPolynomial(coeffs=[11.8677,0.036716,-1.42889e-05,2.52606e-09,-1.69133e-13,63541.3,-24.8084], Tmin=(1340.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C](C=C)C[CH]C(5042)',
    structure = SMILES('[CH2]C=C([CH2])C[CH]C'),
    E0 = (382.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,290.588,291.042],'cm^-1')),
        HinderedRotor(inertia=(0.087591,'amu*angstrom^2'), symmetry=1, barrier=(5.2381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.359784,'amu*angstrom^2'), symmetry=1, barrier=(21.5439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706872,'amu*angstrom^2'), symmetry=1, barrier=(42.3666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.70557e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0379636,'amu*angstrom^2'), symmetry=1, barrier=(96.5445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.454455,0.0685799,-4.4476e-05,1.48416e-08,-2.02838e-12,46086.8,30.9287], Tmin=(100,'K'), Tmax=(1679.59,'K')), NASAPolynomial(coeffs=[15.4741,0.0328104,-1.25312e-05,2.16202e-09,-1.41088e-13,41041.5,-49.3207], Tmin=(1679.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CC(C)C[CH][CH2](5266)',
    structure = SMILES('[CH]=CC(C)C[CH][CH2]'),
    E0 = (550.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,1931.03,4000],'cm^-1')),
        HinderedRotor(inertia=(0.103049,'amu*angstrom^2'), symmetry=1, barrier=(6.78343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386507,'amu*angstrom^2'), symmetry=1, barrier=(25.4422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181731,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784879,'amu*angstrom^2'), symmetry=1, barrier=(51.6668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10305,'amu*angstrom^2'), symmetry=1, barrier=(6.78341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985406,0.0652058,-4.21291e-05,1.45001e-08,-2.13065e-12,66373.4,32.3057], Tmin=(100,'K'), Tmax=(1499.08,'K')), NASAPolynomial(coeffs=[10.7591,0.0391266,-1.60337e-05,2.89497e-09,-1.95272e-13,63443.1,-18.8035], Tmin=(1499.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])CC[CH2](5267)',
    structure = SMILES('[CH]=CC([CH2])CC[CH2]'),
    E0 = (561.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,180,2610.4],'cm^-1')),
        HinderedRotor(inertia=(0.740721,'amu*angstrom^2'), symmetry=1, barrier=(17.0306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740279,'amu*angstrom^2'), symmetry=1, barrier=(17.0205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953337,'amu*angstrom^2'), symmetry=1, barrier=(2.19191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00143749,'amu*angstrom^2'), symmetry=1, barrier=(2.24187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.53058,'amu*angstrom^2'), symmetry=1, barrier=(81.1751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.325951,0.0721471,-5.35778e-05,2.14721e-08,-3.53379e-12,67683.1,33.7321], Tmin=(100,'K'), Tmax=(1428.79,'K')), NASAPolynomial(coeffs=[14.4959,0.0324775,-1.19313e-05,2.04018e-09,-1.33737e-13,63633.9,-39.686], Tmin=(1428.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([C]=C)C[CH]C(5044)',
    structure = SMILES('[CH2]C([C]=C)C[CH]C'),
    E0 = (541.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06387,0.0664489,-4.92513e-05,2.25784e-08,-4.6965e-12,65235.8,32.1006], Tmin=(100,'K'), Tmax=(1078.99,'K')), NASAPolynomial(coeffs=[7.12817,0.043968,-1.79994e-05,3.26958e-09,-2.22774e-13,63927.1,2.38272], Tmin=(1078.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C[CH]C(5045)',
    structure = SMILES('[CH]=CC([CH2])C[CH]C'),
    E0 = (550.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,417.678,1864.07],'cm^-1')),
        HinderedRotor(inertia=(0.00314805,'amu*angstrom^2'), symmetry=1, barrier=(7.76237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627023,'amu*angstrom^2'), symmetry=1, barrier=(7.76237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627024,'amu*angstrom^2'), symmetry=1, barrier=(7.76238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627024,'amu*angstrom^2'), symmetry=1, barrier=(7.76237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56098,'amu*angstrom^2'), symmetry=1, barrier=(69.4476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.877181,0.067174,-4.77127e-05,1.94917e-08,-3.45455e-12,66358.3,32.6932], Tmin=(100,'K'), Tmax=(1278.3,'K')), NASAPolynomial(coeffs=[9.60172,0.0398731,-1.56764e-05,2.78366e-09,-1.86871e-13,64127.8,-11.5396], Tmin=(1278.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(Cds_P)"""),
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
    E0 = (508.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (508.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (668.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (668.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (671.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (912.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (975.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (974.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (963.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (963.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (517.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (572.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (572.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (799.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (849.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (800.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (633.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (578.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (596.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (546.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (562.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (629.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (621.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (595.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (603.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (660.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (508.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (759.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (918.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (799.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (915.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (958.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (967.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (639.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (646.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (626.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (693.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (625.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (653.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (653.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (640.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (595.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (594.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (684.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (683.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['allyl(82)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C(C=C)CC=C(4875)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2][CH]CC[CH]C=C(4160)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2][CH]C[CH]CC=C(5244)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH2])C([CH2])C=C(4122)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH2](502)', '[CH2]C([CH2])C=C(496)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C]CC([CH2])C=C(5245)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C(C=C)C[CH][CH2](5246)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH][CH]CC([CH2])C=C(5247)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C1CC(C=C)C1(5237)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]CCC(=C)C=C(5248)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C=CC(C)C=C(5249)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C[C]([CH2])C[CH2](5250)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH][CH]C([CH2])C[CH2](5251)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH][CH]C([CH2])[CH]C(4150)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2][CH]CC1[CH]CC1(5252)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C1[CH]CC([CH2])C1(5227)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.45491e+07,'s^-1'), n=1.06599, Ea=(69.5416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C1[CH]CC[CH]C1(5253)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2][CH]CC1CC1[CH2](5254)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C1CC([CH2])C1[CH2](5209)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(999564,'s^-1'), n=1.52333, Ea=(53.4157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C1C[CH]CC1[CH2](5255)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][CH]CC(=C)C=C(5256)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C=CC([CH2])C=C(5257)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH2](497)', 'butadiene13(1350)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0273886,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C2H3(60)', '[CH2][CH]CC=C(508)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 432 used for Cds-CsH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CsH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['allyl(82)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.54197,'m^3/(mol*s)'), n=1.98996, Ea=(77.2993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from 72.8 to 77.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.3e+09,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C2H3(60)', '[CH2][CH]C[CH][CH2](510)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][CH]C[C]([CH2])C=C(5258)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][CH][CH]C([CH2])C=C(5259)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH]CC([CH2])[C]=C(5260)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH]=CC([CH2])C[CH][CH2](5261)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]CC(C)=C[CH2](4104)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C[CH]C([CH2])C=C(5262)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2][CH]C(C)C=C[CH2](4105)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH]CC(C)[C]=C(5263)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]CC[C]([CH2])C=C(5264)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2]C([CH][CH]C)C=C(5043)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CCC([CH2])[C]=C(5265)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.572e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]CC([CH2])C=C(4148)'],
    products = ['[CH2][C](C=C)C[CH]C(5042)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=CC(C)C[CH][CH2](5266)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=CC([CH2])CC[CH2](5267)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([C]=C)C[CH]C(5044)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=CC([CH2])C[CH]C(5045)'],
    products = ['[CH2][CH]CC([CH2])C=C(4148)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_4;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1317',
    isomers = [
        '[CH2][CH]CC([CH2])C=C(4148)',
    ],
    reactants = [
        ('allyl(82)', 'butadiene13(1350)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1317',
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


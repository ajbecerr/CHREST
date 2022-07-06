species(
    label = '[CH2]C([CH2])[CH][C]=O(2582)',
    structure = SMILES('[CH2]C([CH2])[CH][C]=O'),
    E0 = (477.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,318.813],'cm^-1')),
        HinderedRotor(inertia=(0.0281889,'amu*angstrom^2'), symmetry=1, barrier=(69.0787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00379562,'amu*angstrom^2'), symmetry=1, barrier=(9.32132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281264,'amu*angstrom^2'), symmetry=1, barrier=(69.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133159,'amu*angstrom^2'), symmetry=1, barrier=(9.3129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41759,0.0538488,-5.61778e-05,3.36165e-08,-7.9844e-12,57553.1,26.3147], Tmin=(100,'K'), Tmax=(1127.64,'K')), NASAPolynomial(coeffs=[9.9417,0.0197888,-5.78525e-06,8.17611e-10,-4.62739e-14,55873.7,-14.7554], Tmin=(1127.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH][CH][C]=O(9637)',
    structure = SMILES('[CH2][CH]C=[C][O]'),
    E0 = (502.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,458.276,458.348],'cm^-1')),
        HinderedRotor(inertia=(0.000816836,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000816843,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47043,0.0263284,1.39799e-06,-2.1241e-08,1.01461e-11,60468.5,21.3394], Tmin=(100,'K'), Tmax=(980.457,'K')), NASAPolynomial(coeffs=[9.92325,0.0131468,-4.7862e-06,8.81255e-10,-6.33441e-14,58179.2,-18.6908], Tmin=(980.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJ) + radical(C=CJO)"""),
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
    label = '[CH2]C([CH2])[C][C]=O(9809)',
    structure = SMILES('[CH2]C([CH2])[C][C]=O'),
    E0 = (758.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,180,1164.67],'cm^-1')),
        HinderedRotor(inertia=(0.0061557,'amu*angstrom^2'), symmetry=1, barrier=(5.93723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257638,'amu*angstrom^2'), symmetry=1, barrier=(5.9236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257948,'amu*angstrom^2'), symmetry=1, barrier=(5.93074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74058,'amu*angstrom^2'), symmetry=1, barrier=(63.0114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33432,0.0616237,-8.70494e-05,6.72732e-08,-2.0086e-11,91309.5,26.7625], Tmin=(100,'K'), Tmax=(958.22,'K')), NASAPolynomial(coeffs=[8.88883,0.0203607,-7.22876e-06,1.14499e-09,-6.90815e-14,90308.3,-7.03094], Tmin=(958.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(758.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C([CH2])[CH][C]=O(9810)',
    structure = SMILES('[CH]C([CH2])[CH][C]=O'),
    E0 = (720.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,194.437,1336.7,1336.77,1336.88],'cm^-1')),
        HinderedRotor(inertia=(2.45642,'amu*angstrom^2'), symmetry=1, barrier=(66.1598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00722298,'amu*angstrom^2'), symmetry=1, barrier=(9.15967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00444482,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00445953,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55717,0.0526161,-5.34782e-05,2.72835e-08,-4.55844e-12,86788.8,25.738], Tmin=(100,'K'), Tmax=(845.473,'K')), NASAPolynomial(coeffs=[10.5679,0.0179253,-6.01765e-06,9.67391e-10,-6.13041e-14,84981.4,-17.8992], Tmin=(845.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(720.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1CC1[C]=O(9811)',
    structure = SMILES('[CH2]C1CC1[C]=O'),
    E0 = (260.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13227,0.0293644,2.24641e-05,-5.21324e-08,2.29551e-11,31430.7,21.9654], Tmin=(100,'K'), Tmax=(939.147,'K')), NASAPolynomial(coeffs=[11.9642,0.0170215,-4.99167e-06,8.41685e-10,-5.99564e-14,28281.6,-31.7845], Tmin=(939.147,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C][CH]C1CC1(9812)',
    structure = SMILES('O=[C][CH]C1CC1'),
    E0 = (230.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15714,0.0274744,2.97277e-05,-6.02329e-08,2.57782e-11,27760.3,20.3478], Tmin=(100,'K'), Tmax=(945.644,'K')), NASAPolynomial(coeffs=[12.4864,0.0165576,-4.94481e-06,8.62255e-10,-6.31069e-14,24341.2,-36.6562], Tmin=(945.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=C)C[C]=O(2579)',
    structure = SMILES('[CH2]C(=C)C[C]=O'),
    E0 = (179.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.051513,'amu*angstrom^2'), symmetry=1, barrier=(18.7985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0518636,'amu*angstrom^2'), symmetry=1, barrier=(18.7912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.050901,'amu*angstrom^2'), symmetry=1, barrier=(18.7872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66359,0.0415019,-1.03508e-05,-1.46014e-08,7.52939e-12,21653.2,21.2212], Tmin=(100,'K'), Tmax=(1094.22,'K')), NASAPolynomial(coeffs=[13.0594,0.0203007,-9.33105e-06,1.86334e-09,-1.36066e-13,17934.6,-40.3792], Tmin=(1094.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(C)=C[C]=O(9670)',
    structure = SMILES('[CH2]C(C)=C[C]=O'),
    E0 = (155.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.835105,'amu*angstrom^2'), symmetry=1, barrier=(19.2007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834119,'amu*angstrom^2'), symmetry=1, barrier=(19.178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.833719,'amu*angstrom^2'), symmetry=1, barrier=(19.1688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54945,0.0590457,-6.39396e-05,3.95624e-08,-1.04434e-11,18810.4,19.6195], Tmin=(100,'K'), Tmax=(894.506,'K')), NASAPolynomial(coeffs=[8.09987,0.029754,-1.48207e-05,2.95469e-09,-2.12236e-13,17638.5,-11.2523], Tmin=(894.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C([CH2])=C[C]=O(9813)',
    structure = SMILES('[CH2]C([CH2])=C[C]=O'),
    E0 = (307.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.12597,'amu*angstrom^2'), symmetry=1, barrier=(25.8882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12599,'amu*angstrom^2'), symmetry=1, barrier=(25.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.126,'amu*angstrom^2'), symmetry=1, barrier=(25.889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50488,0.0572345,-5.84403e-05,3.09075e-08,-6.68244e-12,37035.6,20.0483], Tmin=(100,'K'), Tmax=(1099.74,'K')), NASAPolynomial(coeffs=[11.0466,0.022529,-1.11032e-05,2.21137e-09,-1.58996e-13,34936.9,-26.8919], Tmin=(1099.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=CC=[C][O](9589)',
    structure = SMILES('[CH2]C=C[C]=O'),
    E0 = (194.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.05258,'amu*angstrom^2'), symmetry=1, barrier=(24.201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05258,'amu*angstrom^2'), symmetry=1, barrier=(24.201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33777,0.0398735,-3.70418e-05,1.81131e-08,-3.74289e-12,23480.8,16.0206], Tmin=(100,'K'), Tmax=(1118.07,'K')), NASAPolynomial(coeffs=[7.9573,0.0197688,-1.00692e-05,2.03004e-09,-1.46689e-13,22224.2,-11.7174], Tmin=(1118.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2][C]([CH2])[CH][C]=O(9814)',
    structure = SMILES('[CH2][C]([CH2])[CH][C]=O'),
    E0 = (663.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,1333.91],'cm^-1')),
        HinderedRotor(inertia=(0.00530312,'amu*angstrom^2'), symmetry=1, barrier=(6.69625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291319,'amu*angstrom^2'), symmetry=1, barrier=(6.698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00530571,'amu*angstrom^2'), symmetry=1, barrier=(6.69619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00530272,'amu*angstrom^2'), symmetry=1, barrier=(6.69682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65735,0.0593961,-9.62223e-05,8.63572e-08,-2.88824e-11,79835.8,28.6729], Tmin=(100,'K'), Tmax=(942.726,'K')), NASAPolynomial(coeffs=[4.41772,0.0264986,-1.01696e-05,1.66575e-09,-1.01742e-13,80256.8,20.5114], Tmin=(942.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(663.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(CCJCHO) + radical(Isobutyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]([CH2])C[C]=O(2581)',
    structure = SMILES('[CH2][C]([CH2])C[C]=O'),
    E0 = (495.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,2947.93],'cm^-1')),
        HinderedRotor(inertia=(0.29832,'amu*angstrom^2'), symmetry=1, barrier=(6.85897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298177,'amu*angstrom^2'), symmetry=1, barrier=(6.85568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0116121,'amu*angstrom^2'), symmetry=1, barrier=(71.5411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0116094,'amu*angstrom^2'), symmetry=1, barrier=(71.5286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51959,0.0654283,-0.00011243,1.07044e-07,-3.7672e-11,59688.8,27.3278], Tmin=(100,'K'), Tmax=(916.679,'K')), NASAPolynomial(coeffs=[2.65027,0.0328276,-1.38121e-05,2.39797e-09,-1.53007e-13,60643.9,28.3117], Tmin=(916.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](C)[CH][C]=O(9671)',
    structure = SMILES('[CH2][C](C)[CH][C]=O'),
    E0 = (458.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1589.33],'cm^-1')),
        HinderedRotor(inertia=(0.00245363,'amu*angstrom^2'), symmetry=1, barrier=(4.39197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190413,'amu*angstrom^2'), symmetry=1, barrier=(4.37797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190531,'amu*angstrom^2'), symmetry=1, barrier=(4.38068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190752,'amu*angstrom^2'), symmetry=1, barrier=(4.38576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67412,0.0572992,-8.23064e-05,7.21386e-08,-2.4626e-11,55171,26.1383], Tmin=(100,'K'), Tmax=(900.811,'K')), NASAPolynomial(coeffs=[4.43955,0.0301929,-1.2481e-05,2.19115e-09,-1.42751e-13,55274.3,16.4246], Tmin=(900.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]([CH2])[CH]C=O(9815)',
    structure = SMILES('[CH2]C([CH2])=C[CH][O]'),
    E0 = (423.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71293,0.0410535,-7.93262e-06,-1.72198e-08,8.79803e-12,51074.4,23.7143], Tmin=(100,'K'), Tmax=(1041.05,'K')), NASAPolynomial(coeffs=[11.7807,0.0221058,-9.06762e-06,1.71694e-09,-1.22461e-13,47908.8,-30.3982], Tmin=(1041.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C([CH2])C=C=O(2580)',
    structure = SMILES('[CH2]C([CH2])C=C=O'),
    E0 = (287.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.286408,'amu*angstrom^2'), symmetry=1, barrier=(6.58509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28639,'amu*angstrom^2'), symmetry=1, barrier=(6.58468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627702,'amu*angstrom^2'), symmetry=1, barrier=(14.4321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09136,0.0688467,-0.000100078,8.21214e-08,-2.68731e-11,34653.6,23.2114], Tmin=(100,'K'), Tmax=(841.479,'K')), NASAPolynomial(coeffs=[8.52132,0.0264338,-1.18283e-05,2.18568e-09,-1.47942e-13,33654.4,-9.85906], Tmin=(841.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O) + radical(CJC(C)C=C=O)"""),
)

species(
    label = '[CH2][CH]C[CH][C]=O(2616)',
    structure = SMILES('[CH2][CH]C[CH][C]=O'),
    E0 = (476.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,1652.74,1652.79],'cm^-1')),
        HinderedRotor(inertia=(0.1661,'amu*angstrom^2'), symmetry=1, barrier=(8.24258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333532,'amu*angstrom^2'), symmetry=1, barrier=(64.6535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00425184,'amu*angstrom^2'), symmetry=1, barrier=(8.24363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16632,'amu*angstrom^2'), symmetry=1, barrier=(8.24369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8346,0.0511755,-6.13261e-05,4.91033e-08,-1.65214e-11,57360.1,26.7553], Tmin=(100,'K'), Tmax=(839.036,'K')), NASAPolynomial(coeffs=[5.14519,0.0296716,-1.2654e-05,2.30336e-09,-1.55372e-13,57005.9,12.5647], Tmin=(839.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1C=C([O])C1(9816)',
    structure = SMILES('[CH2]C1[CH]C(=O)C1'),
    E0 = (218.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50574,0.0161091,6.51635e-05,-9.6026e-08,3.80145e-11,26343.7,20.2595], Tmin=(100,'K'), Tmax=(947.072,'K')), NASAPolynomial(coeffs=[12.0266,0.0183054,-5.48228e-06,9.84012e-10,-7.41165e-14,22638.4,-35.1964], Tmin=(947.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C=C[O](9817)',
    structure = SMILES('[CH2]C(=C)C=C[O]'),
    E0 = (140.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2283,0.0430307,1.50481e-05,-6.50533e-08,3.23409e-11,17014.6,20.0949], Tmin=(100,'K'), Tmax=(922.983,'K')), NASAPolynomial(coeffs=[20.6596,0.00529201,8.55132e-07,-2.50556e-10,1.25565e-14,11448.2,-82.8157], Tmin=(922.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C#C[O](9818)',
    structure = SMILES('[CH2]C([CH2])[C]=C=O'),
    E0 = (493.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.362801,'amu*angstrom^2'), symmetry=1, barrier=(8.34151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.362475,'amu*angstrom^2'), symmetry=1, barrier=(8.33402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.50701,'amu*angstrom^2'), symmetry=1, barrier=(103.625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02895,0.0729765,-0.000120849,1.06643e-07,-3.61612e-11,59404.3,23.265], Tmin=(100,'K'), Tmax=(872.549,'K')), NASAPolynomial(coeffs=[7.77751,0.0258499,-1.20021e-05,2.21523e-09,-1.48194e-13,58842.9,-4.84133], Tmin=(872.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O) + radical(CJC(C)C=C=O) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = '[CH2]C([CH2])[C]=C[O](9819)',
    structure = SMILES('[CH2]C([CH2])[C]=C[O]'),
    E0 = (532.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,364.844,3276.1],'cm^-1')),
        HinderedRotor(inertia=(0.211785,'amu*angstrom^2'), symmetry=1, barrier=(20.0155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852719,'amu*angstrom^2'), symmetry=1, barrier=(80.4583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.851971,'amu*angstrom^2'), symmetry=1, barrier=(80.4577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.566681,0.0586772,-5.57203e-05,2.69578e-08,-4.8869e-12,64171.3,29.0205], Tmin=(100,'K'), Tmax=(1597.2,'K')), NASAPolynomial(coeffs=[14.8082,0.0119882,-1.52087e-06,1.42371e-11,6.7165e-15,61028,-41.9541], Tmin=(1597.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)[C][C]=O(9666)',
    structure = SMILES('[CH2]C(C)[C][C]=O'),
    E0 = (553.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.110357,'amu*angstrom^2'), symmetry=1, barrier=(2.53733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108669,'amu*angstrom^2'), symmetry=1, barrier=(2.49851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109533,'amu*angstrom^2'), symmetry=1, barrier=(2.51839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.84626,'amu*angstrom^2'), symmetry=1, barrier=(88.4331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54162,0.0569827,-6.26091e-05,3.64919e-08,-7.1099e-12,66636.5,24.9486], Tmin=(100,'K'), Tmax=(707.203,'K')), NASAPolynomial(coeffs=[8.57968,0.0246536,-9.90087e-06,1.75835e-09,-1.1756e-13,65454.1,-7.88998], Tmin=(707.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH2])[C]=[C]O(9820)',
    structure = SMILES('[CH2]C([CH2])[C]=[C]O'),
    E0 = (630.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,465.235],'cm^-1')),
        HinderedRotor(inertia=(0.0504699,'amu*angstrom^2'), symmetry=1, barrier=(7.75185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0504698,'amu*angstrom^2'), symmetry=1, barrier=(7.75183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0504704,'amu*angstrom^2'), symmetry=1, barrier=(7.75185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.508935,'amu*angstrom^2'), symmetry=1, barrier=(78.1688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790012,0.0626319,-7.03485e-05,4.05697e-08,-8.82004e-12,75976.2,29.9213], Tmin=(100,'K'), Tmax=(1298.43,'K')), NASAPolynomial(coeffs=[14.2743,0.0121495,-1.69882e-06,1.80813e-11,9.01322e-15,73228.4,-35.7514], Tmin=(1298.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(630.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]([CH2])C=[C]O(9821)',
    structure = SMILES('[CH2]C(=C)[CH][C]O'),
    E0 = (514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,345.896,345.917,345.933],'cm^-1')),
        HinderedRotor(inertia=(0.322818,'amu*angstrom^2'), symmetry=1, barrier=(27.4121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322758,'amu*angstrom^2'), symmetry=1, barrier=(27.4128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.32278,'amu*angstrom^2'), symmetry=1, barrier=(27.4122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322829,'amu*angstrom^2'), symmetry=1, barrier=(27.4128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615381,0.0677885,-7.03812e-05,3.60459e-08,-7.18049e-12,61947.4,21.9004], Tmin=(100,'K'), Tmax=(1231.26,'K')), NASAPolynomial(coeffs=[16.9548,0.0147062,-5.71254e-06,1.0308e-09,-7.08533e-14,57923.8,-60.3272], Tmin=(1231.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(Allyl_P) + radical(CH2_triplet)"""),
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
    E0 = (477.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (965.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (940.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (970.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (932.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (483.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (485.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (500.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (500.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (526.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (616.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (581.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (874.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (653.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (619.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (630.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (477.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (635.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (486.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (541.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (719.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (672.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (736.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (705.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (781.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (597.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['HCCO(2227)', 'allyl(82)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH][CH2](497)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH][C]=O(9637)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C([CH2])[C][C]=O(9809)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]C([CH2])[CH][C]=O(9810)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C1CC1[C]=O(9811)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.18842e+14,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['O=[C][CH]C1CC1(9812)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C(=C)C[C]=O(2579)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C(C)=C[C]=O(9670)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2(T)(20)', 'C=CC=[C][O](9589)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C][O](6861)', 'allyl(82)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2][C]([CH2])[CH][C]=O(9814)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]([CH2])C[C]=O(2581)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2][C](C)[CH][C]=O(9671)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2][C]([CH2])[CH]C=O(9815)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C([CH2])C=C=O(2580)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2][CH]C[CH][C]=O(2616)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C1C=C([O])C1(9816)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C(=C)C=C[O](9817)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2]C([CH2])C#C[O](9818)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCCO(2227)', '[CH2][CH][CH2](497)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH2])[C]=C[O](9819)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C)[C][C]=O(9666)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])[C]=[C]O(9820)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]([CH2])C=[C]O(9821)'],
    products = ['[CH2]C([CH2])[CH][C]=O(2582)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2332',
    isomers = [
        '[CH2]C([CH2])[CH][C]=O(2582)',
    ],
    reactants = [
        ('HCCO(2227)', 'allyl(82)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2332',
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


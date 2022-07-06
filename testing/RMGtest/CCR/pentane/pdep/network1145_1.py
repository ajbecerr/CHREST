species(
    label = 'C[CH][CH]CC(C)[CH]C(1248)',
    structure = SMILES('C[CH][CH]CC(C)[CH]C'),
    E0 = (341.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,211.095,896.703,2754.69,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79078,0.0691751,-3.0975e-05,5.72233e-09,-3.69667e-13,41152.1,33.9842], Tmin=(100,'K'), Tmax=(2677.54,'K')), NASAPolynomial(coeffs=[41.8126,0.0185235,-7.71796e-06,1.20621e-09,-6.69998e-14,16444.7,-204.632], Tmin=(2677.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
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
    label = 'butene2t(396)',
    structure = SMILES('CC=CC'),
    E0 = (-28.2774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.464594,'amu*angstrom^2'), symmetry=1, barrier=(10.6819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460586,'amu*angstrom^2'), symmetry=1, barrier=(10.5898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52,0.0276811,9.44413e-07,-1.25016e-08,4.67173e-12,-3343.75,13.0195], Tmin=(100,'K'), Tmax=(1157.77,'K')), NASAPolynomial(coeffs=[6.04212,0.0260075,-1.04845e-05,1.90884e-09,-1.30603e-13,-4862.7,-7.52639], Tmin=(1157.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene2t""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH]C(C)CC=CC(1244)',
    structure = SMILES('C[CH]C(C)CC=CC'),
    E0 = (64.9597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20707,0.0812628,-4.58987e-05,1.24449e-08,-1.34478e-12,7974.21,36.0539], Tmin=(100,'K'), Tmax=(2106.07,'K')), NASAPolynomial(coeffs=[23.2938,0.0366295,-1.41106e-05,2.38285e-09,-1.50403e-13,-1924.98,-94.8291], Tmin=(2106.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.9597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S)"""),
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
    label = 'C[CH][CH]CC[CH]C(328)',
    structure = SMILES('C[CH][CH]CC[CH]C'),
    E0 = (367.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,875.63,1946.17,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3350.7,'J/mol'), sigma=(6.3658,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=523.37 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09928,0.0728686,-8.23053e-05,7.88527e-08,-3.17483e-11,44320.4,34.783], Tmin=(100,'K'), Tmax=(826.938,'K')), NASAPolynomial(coeffs=[-1.37392,0.0640962,-2.87802e-05,5.37859e-09,-3.68367e-13,45438.5,50.5315], Tmin=(826.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]C[CH]C(C)C(4492)',
    structure = SMILES('C[CH][CH]C[CH]C(C)C'),
    E0 = (338.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79078,0.0691751,-3.0975e-05,5.72233e-09,-3.69667e-13,40749.5,33.291], Tmin=(100,'K'), Tmax=(2677.54,'K')), NASAPolynomial(coeffs=[41.8126,0.0185235,-7.71796e-06,1.20621e-09,-6.69998e-14,16042.2,-205.325], Tmin=(2677.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]C)C(C)[CH]C(1230)',
    structure = SMILES('[CH2]C([CH]C)C(C)[CH]C'),
    E0 = (350.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,210.057,640.305,3082.52],'cm^-1')),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779289,'amu*angstrom^2'), symmetry=1, barrier=(2.32772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.237081,0.0808306,-4.72681e-05,1.4049e-08,-1.74851e-12,42237,38.4149], Tmin=(100,'K'), Tmax=(1737.37,'K')), NASAPolynomial(coeffs=[14.0427,0.0490453,-1.98253e-05,3.51852e-09,-2.33214e-13,37439.9,-35.8151], Tmin=(1737.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'CHCH3(T)(21)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82365,-0.000909765,3.21389e-05,-3.73492e-08,1.33096e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.802,'K')), NASAPolynomial(coeffs=[4.3048,0.00943081,-3.27566e-06,5.95138e-10,-4.27321e-14,40709.2,1.84242], Tmin=(960.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH]C[CH]C(160)',
    structure = SMILES('C[CH][CH]C[CH]C'),
    E0 = (391.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,535.998,2986.3,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0181257,'amu*angstrom^2'), symmetry=1, barrier=(3.33465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0181257,'amu*angstrom^2'), symmetry=1, barrier=(3.33465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0181257,'amu*angstrom^2'), symmetry=1, barrier=(3.33465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0181257,'amu*angstrom^2'), symmetry=1, barrier=(3.33465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0181257,'amu*angstrom^2'), symmetry=1, barrier=(3.33465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77602,0.0577092,-6.69897e-05,6.85215e-08,-2.84115e-11,47156.5,30.1092], Tmin=(100,'K'), Tmax=(840.292,'K')), NASAPolynomial(coeffs=[-2.35256,0.0557803,-2.50208e-05,4.65912e-09,-3.17789e-13,48612.2,53.8426], Tmin=(840.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH][CH]C(3874)',
    structure = SMILES('[CH][CH]C'),
    E0 = (522.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1852.85,1853.13,1853.16],'cm^-1')),
        HinderedRotor(inertia=(0.0369856,'amu*angstrom^2'), symmetry=1, barrier=(8.28142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370949,'amu*angstrom^2'), symmetry=1, barrier=(8.27684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19279,0.0163284,-1.76328e-06,-3.55309e-09,1.20622e-12,62877.3,14.3587], Tmin=(100,'K'), Tmax=(1426.65,'K')), NASAPolynomial(coeffs=[5.59458,0.0149965,-6.04303e-06,1.1011e-09,-7.44871e-14,61642.3,-0.00873183], Tmin=(1426.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)[CH]C(114)',
    structure = SMILES('[CH2]C(C)[CH]C'),
    E0 = (225.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1290.57],'cm^-1')),
        HinderedRotor(inertia=(0.100358,'amu*angstrom^2'), symmetry=1, barrier=(2.30744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10319,'amu*angstrom^2'), symmetry=1, barrier=(2.37253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0987185,'amu*angstrom^2'), symmetry=1, barrier=(2.26973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971084,'amu*angstrom^2'), symmetry=1, barrier=(2.23271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80047,0.0421184,-1.07295e-05,-7.76584e-09,4.08549e-12,27247.4,23.8047], Tmin=(100,'K'), Tmax=(1119,'K')), NASAPolynomial(coeffs=[7.75174,0.0320963,-1.23774e-05,2.20164e-09,-1.48914e-13,25211,-8.72346], Tmin=(1119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH]C(C)C[CH][CH]C(3930)',
    structure = SMILES('[CH]C(C)C[CH][CH]C'),
    E0 = (615.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,219.479,1133.55,1226.17,1412.49,1676.71,1962.91],'cm^-1')),
        HinderedRotor(inertia=(0.0895117,'amu*angstrom^2'), symmetry=1, barrier=(2.85145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895117,'amu*angstrom^2'), symmetry=1, barrier=(2.85145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895117,'amu*angstrom^2'), symmetry=1, barrier=(2.85145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895117,'amu*angstrom^2'), symmetry=1, barrier=(2.85145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895117,'amu*angstrom^2'), symmetry=1, barrier=(2.85145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895117,'amu*angstrom^2'), symmetry=1, barrier=(2.85145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955494,0.071958,-6.73597e-05,5.02407e-08,-1.79415e-11,74166.2,33.2951], Tmin=(100,'K'), Tmax=(755.539,'K')), NASAPolynomial(coeffs=[3.43643,0.0545043,-2.41335e-05,4.53298e-09,-3.13752e-13,73914.6,22.8372], Tmin=(755.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC(C)[CH]C(1269)',
    structure = SMILES('[CH][CH]CC(C)[CH]C'),
    E0 = (619.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,211.398,976.848,1097.8,1324.35,1554.74,1953.7],'cm^-1')),
        HinderedRotor(inertia=(0.124017,'amu*angstrom^2'), symmetry=1, barrier=(3.30786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124017,'amu*angstrom^2'), symmetry=1, barrier=(3.30786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124017,'amu*angstrom^2'), symmetry=1, barrier=(3.30786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124017,'amu*angstrom^2'), symmetry=1, barrier=(3.30786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124017,'amu*angstrom^2'), symmetry=1, barrier=(3.30786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124017,'amu*angstrom^2'), symmetry=1, barrier=(3.30786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26337,0.0660642,-4.05541e-05,1.35985e-08,-2.08965e-12,74567.1,32.1677], Tmin=(100,'K'), Tmax=(1317.39,'K')), NASAPolynomial(coeffs=[6.71103,0.0495234,-2.17204e-05,4.0677e-09,-2.80993e-13,73131.7,4.38421], Tmin=(1317.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    label = 'C[C]C(C)C[CH][CH]C(4493)',
    structure = SMILES('C[C]C(C)C[CH][CH]C'),
    E0 = (595.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40278,0.0689683,1.46379e-05,-1.49281e-07,1.36225e-10,71685,33.0947], Tmin=(100,'K'), Tmax=(449.471,'K')), NASAPolynomial(coeffs=[4.08657,0.063459,-2.82971e-05,5.35449e-09,-3.73219e-13,71258.1,20.2284], Tmin=(449.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH][C]CC(C)[CH]C(4494)',
    structure = SMILES('C[CH][C]CC(C)[CH]C'),
    E0 = (595.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619359,0.0807149,-5.35665e-05,2.02359e-08,-3.53191e-12,71730,35.6346], Tmin=(100,'K'), Tmax=(1202.57,'K')), NASAPolynomial(coeffs=[7.2972,0.0585031,-2.58613e-05,4.87705e-09,-3.39012e-13,70123.9,2.18599], Tmin=(1202.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C][CH]CC(C)[CH]C(4495)',
    structure = SMILES('C[C][CH]CC(C)[CH]C'),
    E0 = (595.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620547,0.0813101,-5.79089e-05,2.60453e-08,-5.75828e-12,71731.1,35.4806], Tmin=(100,'K'), Tmax=(960.052,'K')), NASAPolynomial(coeffs=[5.19167,0.0622653,-2.81541e-05,5.38398e-09,-3.78185e-13,70853.4,13.6138], Tmin=(960.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]C1CC(C)C1C(4487)',
    structure = SMILES('C[CH]C1CC(C)C1C'),
    E0 = (86.1697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.697114,0.0506046,5.68393e-05,-1.00578e-07,3.93198e-11,10502.2,30.3364], Tmin=(100,'K'), Tmax=(997.651,'K')), NASAPolynomial(coeffs=[15.027,0.0468005,-1.81064e-05,3.40664e-09,-2.45259e-13,4972.98,-52.1447], Tmin=(997.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.1697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]CCC(C)=CC(4496)',
    structure = SMILES('C[CH]CCC(C)=CC'),
    E0 = (55.2834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687643,0.0787071,-4.45456e-05,1.28004e-08,-1.58469e-12,6761.85,31.7221], Tmin=(100,'K'), Tmax=(1639.54,'K')), NASAPolynomial(coeffs=[10.0282,0.0559188,-2.36968e-05,4.32291e-09,-2.92028e-13,3699,-17.9589], Tmin=(1639.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.2834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C=CC(C)CC(4031)',
    structure = SMILES('C[CH]C=CC(C)CC'),
    E0 = (10.3874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0561617,0.0717752,-7.84247e-06,-3.22197e-08,1.50391e-11,1404.11,32.0931], Tmin=(100,'K'), Tmax=(1052.76,'K')), NASAPolynomial(coeffs=[14.6734,0.0479564,-1.91004e-05,3.52995e-09,-2.46965e-13,-3431.34,-47.5265], Tmin=(1052.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.3874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S)"""),
)

species(
    label = 'C=CC(C)CC[CH]C(3801)',
    structure = SMILES('C=CC(C)CC[CH]C'),
    E0 = (74.8316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.17103,0.0816069,-4.70021e-05,1.34621e-08,-1.56758e-12,9158.72,36.5529], Tmin=(100,'K'), Tmax=(1929.68,'K')), NASAPolynomial(coeffs=[19.1699,0.0415151,-1.58372e-05,2.69515e-09,-1.72662e-13,1694.41,-69.4697], Tmin=(1929.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]CC(C)=CC(4497)',
    structure = SMILES('C[CH][CH]CC(C)=CC'),
    E0 = (249.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,180,400,3152.31],'cm^-1')),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17031,0.0665187,-3.05677e-05,6.00427e-09,-4.422e-13,30073.5,28.8102], Tmin=(100,'K'), Tmax=(3278.27,'K')), NASAPolynomial(coeffs=[48.8981,0.0095073,-4.48342e-06,7.00143e-10,-3.77353e-14,-565.886,-252.108], Tmin=(3278.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C=CC(C)[CH]C(4498)',
    structure = SMILES('C[CH]C=CC(C)[CH]C'),
    E0 = (204.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215087,0.0708218,-1.96097e-05,-1.43963e-08,7.68642e-12,24793.9,34.0858], Tmin=(100,'K'), Tmax=(1136.31,'K')), NASAPolynomial(coeffs=[13.4909,0.0474171,-1.95089e-05,3.61161e-09,-2.50456e-13,20270.7,-38.2858], Tmin=(1136.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'C=CC(C)C[CH][CH]C(3837)',
    structure = SMILES('C=CC(C)C[CH][CH]C'),
    E0 = (269.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,215.829,823.733,1222.4,1702.25],'cm^-1')),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15014,0.0706793,-3.54891e-05,8.33376e-09,-7.87455e-13,32481.1,34.2808], Tmin=(100,'K'), Tmax=(2185.09,'K')), NASAPolynomial(coeffs=[14.1606,0.0468625,-1.91396e-05,3.34557e-09,-2.16748e-13,26795.3,-38.6569], Tmin=(2185.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C=C[CH]CC(C)[CH]C(4499)',
    structure = SMILES('[CH2]C=CCC(C)[CH]C'),
    E0 = (216.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0834083,0.077896,-3.66884e-05,1.0662e-09,2.82546e-12,26190.8,35.8462], Tmin=(100,'K'), Tmax=(1202.06,'K')), NASAPolynomial(coeffs=[15.0069,0.0455493,-1.86212e-05,3.41191e-09,-2.34185e-13,21272,-45.1024], Tmin=(1202.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]C(3856)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (450.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2810.32],'cm^-1')),
        HinderedRotor(inertia=(0.00169393,'amu*angstrom^2'), symmetry=1, barrier=(9.48541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010072,'amu*angstrom^2'), symmetry=1, barrier=(56.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100766,'amu*angstrom^2'), symmetry=1, barrier=(56.4312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79546,0.0335488,-4.51222e-05,5.0661e-08,-2.13657e-11,54172.4,20.3345], Tmin=(100,'K'), Tmax=(854.972,'K')), NASAPolynomial(coeffs=[-1.26875,0.0345785,-1.53757e-05,2.86235e-09,-1.94746e-13,55524.7,43.1492], Tmin=(854.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]CC=CC(3959)',
    structure = SMILES('C[CH][CH]CC=CC'),
    E0 = (288.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,199.254,410.84,3100.51],'cm^-1')),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21183,0.0544927,-2.28084e-05,3.60194e-09,-1.42453e-13,34783.4,28.0308], Tmin=(100,'K'), Tmax=(2468.19,'K')), NASAPolynomial(coeffs=[28.511,0.0216028,-8.73419e-06,1.39784e-09,-8.1001e-14,18837,-128.613], Tmin=(2468.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]C(1186)',
    structure = SMILES('C[CH][CH]C'),
    E0 = (244.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2670.04],'cm^-1')),
        HinderedRotor(inertia=(0.341695,'amu*angstrom^2'), symmetry=1, barrier=(7.85623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155291,'amu*angstrom^2'), symmetry=1, barrier=(7.85531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136081,'amu*angstrom^2'), symmetry=1, barrier=(68.8251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62536,0.022997,-1.74547e-06,-3.25078e-09,7.61424e-13,29450.7,14.8448], Tmin=(100,'K'), Tmax=(2137.09,'K')), NASAPolynomial(coeffs=[11.4859,0.0201914,-8.13358e-06,1.34908e-09,-8.16544e-14,23371.9,-35.4092], Tmin=(2137.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]C[CH][CH]C(3962)',
    structure = SMILES('C[CH][CH]C[CH][CH]C'),
    E0 = (562.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3527.57,3926.39,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33443,0.0735713,-0.000107041,1.15927e-07,-4.73572e-11,67694.3,35.8086], Tmin=(100,'K'), Tmax=(858.862,'K')), NASAPolynomial(coeffs=[-5.5633,0.0679124,-3.11683e-05,5.81031e-09,-3.9425e-13,70272.7,76.1494], Tmin=(858.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]C[C](C)[CH]C(4500)',
    structure = SMILES('C[CH][CH]C[C](C)[CH]C'),
    E0 = (527.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1765.94,2581.8,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.4931,0.091778,-0.000127515,1.30183e-07,-5.22604e-11,63507.2,39.8762], Tmin=(100,'K'), Tmax=(837.616,'K')), NASAPolynomial(coeffs=[-3.40655,0.0753008,-3.51508e-05,6.64087e-09,-4.55645e-13,65391.8,65.3488], Tmin=(837.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH][CH]C(C)[CH]C(4501)',
    structure = SMILES('C[CH][CH][CH]C(C)[CH]C'),
    E0 = (536.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,1698.15,2965.81,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658945,0.0823981,-8.94173e-05,8.31501e-08,-3.39568e-11,64603.6,40.7044], Tmin=(100,'K'), Tmax=(790.11,'K')), NASAPolynomial(coeffs=[-0.589045,0.0706006,-3.26283e-05,6.21492e-09,-4.31746e-13,65366.3,50.0095], Tmin=(790.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH][CH]C(3834)',
    structure = SMILES('[CH2]C([CH]C)C[CH][CH]C'),
    E0 = (546.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180,1264.05,3750.06,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442845,0.0857364,-9.15222e-05,7.75899e-08,-2.89865e-11,65880.5,40.7089], Tmin=(100,'K'), Tmax=(812.499,'K')), NASAPolynomial(coeffs=[1.77095,0.0654815,-2.88057e-05,5.35245e-09,-3.6637e-13,66117.4,37.3635], Tmin=(812.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)C[CH][CH]C(3839)',
    structure = SMILES('[CH2][CH]C(C)C[CH][CH]C'),
    E0 = (546.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,411.584,1157.99,3528.67,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)[CH]C(4502)',
    structure = SMILES('[CH2][CH][CH]CC(C)[CH]C'),
    E0 = (546.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,411.584,1157.99,3528.67,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[C](C)CC(4039)',
    structure = SMILES('C[CH][CH]C[C](C)CC'),
    E0 = (332.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98806,0.0642371,-2.55546e-05,3.43802e-09,-6.58431e-14,39989.8,28.123], Tmin=(100,'K'), Tmax=(2720.54,'K')), NASAPolynomial(coeffs=[54.396,0.00666161,-3.73934e-06,5.25479e-10,-2.18015e-14,5353.71,-283.59], Tmin=(2720.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH]C(C)[CH]C(4503)',
    structure = SMILES('C[CH]C[CH]C(C)[CH]C'),
    E0 = (341.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,461.027,1409.86,1595.87,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524493,'amu*angstrom^2'), symmetry=1, barrier=(5.50824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20779,0.0730421,-3.5651e-05,7.64594e-09,-6.26054e-13,41191.7,36.1129], Tmin=(100,'K'), Tmax=(2816.47,'K')), NASAPolynomial(coeffs=[33.616,0.0270171,-1.114e-05,1.84432e-09,-1.11102e-13,22935.6,-153.797], Tmin=(2816.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC(C)C[CH][CH]C(4009)',
    structure = SMILES('[CH2]CC(C)C[CH][CH]C'),
    E0 = (352.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.43,716.6,2010.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91664,0.0586225,7.61378e-05,-2.97843e-07,2.71085e-10,42442.7,33.283], Tmin=(100,'K'), Tmax=(393.52,'K')), NASAPolynomial(coeffs=[2.95224,0.0672715,-2.99223e-05,5.66129e-09,-3.94668e-13,42212.7,27.3662], Tmin=(393.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]CC(C)[CH]C(1249)',
    structure = SMILES('[CH2]C[CH]CC(C)[CH]C'),
    E0 = (352.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,215.49,821.979,1188.01,3212.6],'cm^-1')),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18682,0.0746523,-3.81666e-05,8.96888e-09,-8.44223e-13,42479.7,35.2217], Tmin=(100,'K'), Tmax=(2166.95,'K')), NASAPolynomial(coeffs=[14.2593,0.0505217,-2.14631e-05,3.83006e-09,-2.51362e-13,36814.1,-37.9549], Tmin=(2166.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C(C)CC(4041)',
    structure = SMILES('C[CH][CH][CH]C(C)CC'),
    E0 = (341.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,211.095,896.703,2754.69,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860585,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79078,0.0691751,-3.0975e-05,5.72233e-09,-3.69666e-13,41152.1,33.9842], Tmin=(100,'K'), Tmax=(2677.53,'K')), NASAPolynomial(coeffs=[41.8122,0.0185239,-7.71817e-06,1.20624e-09,-6.70023e-14,16445,-204.629], Tmin=(2677.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C(CC)C[CH][CH]C(365)',
    structure = SMILES('[CH2]C(CC)C[CH][CH]C'),
    E0 = (352.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,215.647,509.563,1150.94,3517.27],'cm^-1')),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72881,0.0625006,5.5945e-05,-2.48776e-07,2.23958e-10,42430,33.9051], Tmin=(100,'K'), Tmax=(416.553,'K')), NASAPolynomial(coeffs=[3.63613,0.064977,-2.78436e-05,5.14836e-09,-3.53243e-13,42090.7,24.2084], Tmin=(416.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]CC[C](C)[CH]C(4504)',
    structure = SMILES('C[CH]CC[C](C)[CH]C'),
    E0 = (332.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39918,0.0681473,-3.03074e-05,5.40738e-09,-3.30723e-13,40029.9,30.2754], Tmin=(100,'K'), Tmax=(2841.86,'K')), NASAPolynomial(coeffs=[54.964,0.00554769,-3.27609e-06,4.76096e-10,-2.09482e-14,5555.38,-286.311], Tmin=(2841.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C(C)[CH][CH]CC(1247)',
    structure = SMILES('C[CH]C(C)[CH][CH]CC'),
    E0 = (341.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,406.653,474.659,2309.48,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0393842,'amu*angstrom^2'), symmetry=1, barrier=(3.68198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41927,0.0714365,-3.35488e-05,6.62165e-09,-4.71819e-13,41183.4,35.1861], Tmin=(100,'K'), Tmax=(2585.07,'K')), NASAPolynomial(coeffs=[38.204,0.0221299,-9.35534e-06,1.52148e-09,-8.87472e-14,19621.7,-182.135], Tmin=(2585.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CCC(C)[CH]C(4505)',
    structure = SMILES('[CH2][CH]CCC(C)[CH]C'),
    E0 = (352.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,218.038,1006.99,1394.93,1979.65],'cm^-1')),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866075,0.0771162,-4.19528e-05,1.11256e-08,-1.23983e-12,42495.2,36.5804], Tmin=(100,'K'), Tmax=(1800.21,'K')), NASAPolynomial(coeffs=[10.8904,0.0548426,-2.33937e-05,4.25269e-09,-2.85371e-13,38886,-17.6743], Tmin=(1800.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)CC[CH]C(4044)',
    structure = SMILES('[CH2]C([CH]C)CC[CH]C'),
    E0 = (352.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.219,291.627,1201.19,3327.11],'cm^-1')),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557675,0.080864,-5.17249e-05,1.96292e-08,-3.53158e-12,42491.2,37.735], Tmin=(100,'K'), Tmax=(1156.65,'K')), NASAPolynomial(coeffs=[6.28324,0.0610634,-2.60465e-05,4.8287e-09,-3.32569e-13,41166.7,9.2792], Tmin=(1156.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH][C](C)C[CH]CC(1246)',
    structure = SMILES('C[CH][C](C)C[CH]CC'),
    E0 = (332.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,193.318,206.336,2354.63,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0180474,'amu*angstrom^2'), symmetry=1, barrier=(6.95271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61391,0.0665187,-2.81663e-05,4.36151e-09,-1.72822e-13,40021.3,29.3355], Tmin=(100,'K'), Tmax=(2688.96,'K')), NASAPolynomial(coeffs=[51.0863,0.00987652,-5.19538e-06,8.05051e-10,-4.10047e-14,8362.76,-262.858], Tmin=(2688.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S)"""),
)

species(
    label = '[CH2][CH]C(C)CC[CH]C(4004)',
    structure = SMILES('[CH2][CH]C(C)CC[CH]C'),
    E0 = (352.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,218.038,1006.99,1394.93,1979.65],'cm^-1')),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103829,'amu*angstrom^2'), symmetry=1, barrier=(3.23432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866075,0.0771162,-4.19528e-05,1.11256e-08,-1.23983e-12,42495.2,36.5804], Tmin=(100,'K'), Tmax=(1800.21,'K')), NASAPolynomial(coeffs=[10.8904,0.0548426,-2.33937e-05,4.25269e-09,-2.85371e-13,38886,-17.6743], Tmin=(1800.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH]CC(363)',
    structure = SMILES('[CH2]C([CH]C)C[CH]CC'),
    E0 = (352.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,209.12,746.778,915.592,3452.43],'cm^-1')),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542927,0.0816895,-5.70204e-05,2.68127e-08,-6.3819e-12,42492.8,37.6346], Tmin=(100,'K'), Tmax=(894.836,'K')), NASAPolynomial(coeffs=[4.39928,0.0644514,-2.81247e-05,5.28517e-09,-3.67574e-13,41802.6,19.4584], Tmin=(894.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)C[CH]CC(1250)',
    structure = SMILES('[CH2][CH]C(C)C[CH]CC'),
    E0 = (352.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,215.49,821.979,1188.01,3212.6],'cm^-1')),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718495,'amu*angstrom^2'), symmetry=1, barrier=(2.44252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18682,0.0746523,-3.81666e-05,8.96888e-09,-8.44223e-13,42479.7,35.2217], Tmin=(100,'K'), Tmax=(2166.95,'K')), NASAPolynomial(coeffs=[14.2593,0.0505217,-2.14631e-05,3.83006e-09,-2.51362e-13,36814.1,-37.9549], Tmin=(2166.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)CC(4046)',
    structure = SMILES('[CH2][CH][CH]CC(C)CC'),
    E0 = (352.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.43,716.6,2010.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91664,0.0586225,7.61378e-05,-2.97843e-07,2.71085e-10,42442.7,33.283], Tmin=(100,'K'), Tmax=(393.52,'K')), NASAPolynomial(coeffs=[2.95224,0.0672715,-2.99223e-05,5.66129e-09,-3.94668e-13,42212.7,27.3662], Tmin=(393.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (341.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (341.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (786.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (536.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (507.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (790.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (803.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (807.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (811.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (807.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (807.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (807.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (349.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (405.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (405.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (366.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (470.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (424.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (483.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (436.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (438.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (453.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (403.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (694.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (698.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (738.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (748.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (758.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (758.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (758.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (496.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (478.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (458.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (493.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (464.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (499.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (457.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (488.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (493.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (428.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (459.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (426.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (398.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (395.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (395.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['m1_allyl(186)', 'butene2t(396)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH]C(C)CC=CC(1244)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', 'C[CH][CH]CC[CH]C(328)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH][CH]C[CH]C(C)C(4492)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C([CH]C)C(C)[CH]C(1230)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHCH3(T)(21)', 'C[CH][CH]C[CH]C(160)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH]C(3874)', '[CH2]C(C)[CH]C(114)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH]C(C)C[CH][CH]C(3930)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH][CH]CC(C)[CH]C(1269)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C[C]C(C)C[CH][CH]C(4493)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'C[CH][C]CC(C)[CH]C(4494)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', 'C[C][CH]CC(C)[CH]C(4495)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH]C1CC(C)C1C(4487)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH]CCC(C)=CC(4496)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH]C=CC(C)CC(4031)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C=CC(C)CC[CH]C(3801)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', 'C[CH][CH]CC(C)=CC(4497)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', 'C[CH]C=CC(C)[CH]C(4498)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'C=CC(C)C[CH][CH]C(3837)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', 'C=C[CH]CC(C)[CH]C(4499)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][CH]C(3856)', 'butene2t(396)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00858789,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3(17)', 'C[CH][CH]CC=CC(3959)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(10100,'cm^3/(mol*s)'), n=2.41, Ea=(28.4512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 435 used for Cds-CsH_Cds-CsH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['m1_allyl(186)', 'C[CH][CH]C(1186)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][CH]C(3856)', 'C[CH][CH]C(1186)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
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
    reactants = ['CH3(17)', 'C[CH][CH]C[CH][CH]C(3962)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'C[CH][CH]C[C](C)[CH]C(4500)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'C[CH][CH][CH]C(C)[CH]C(4501)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2][CH]C(C)C[CH][CH]C(3839)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH][CH]C[C](C)CC(4039)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[CH]C[CH]C(C)[CH]C(4503)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['[CH2]CC(C)C[CH][CH]C(4009)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C[CH]CC(C)[CH]C(1249)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH][CH][CH]C(C)CC(4041)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.3913e+06,'s^-1'), n=1.66106, Ea=(123.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH]CC[C](C)[CH]C(4504)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C[CH]C(C)[CH][CH]CC(1247)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]CCC(C)[CH]C(4505)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH]C)CC[CH]C(4044)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C[CH][CH]CC(C)[CH]C(1248)'],
    products = ['C[CH][C](C)C[CH]CC(1246)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C(C)CC[CH]C(4004)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C([CH]C)C[CH]CC(363)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH]C(C)C[CH]CC(1250)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH][CH]CC(C)CC(4046)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1145',
    isomers = [
        'C[CH][CH]CC(C)[CH]C(1248)',
    ],
    reactants = [
        ('m1_allyl(186)', 'butene2t(396)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1145',
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


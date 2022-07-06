species(
    label = '[CH2]C([CH2])C[CH][CH]C(607)',
    structure = SMILES('[CH2]C([CH2])C[CH][CH]C'),
    E0 = (577.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,375.685,1769.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.876276,0.0742873,-7.9907e-05,6.43255e-08,-2.21286e-11,69592,35.4413], Tmin=(100,'K'), Tmax=(870.302,'K')), NASAPolynomial(coeffs=[3.13779,0.05207,-2.12369e-05,3.77352e-09,-2.50661e-13,69646.1,27.4173], Tmin=(870.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH2])CC=CC(604)',
    structure = SMILES('[CH2]C([CH2])CC=CC'),
    E0 = (301.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,276.958,2664.34],'cm^-1')),
        HinderedRotor(inertia=(0.468433,'amu*angstrom^2'), symmetry=1, barrier=(10.7702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.39484,'amu*angstrom^2'), symmetry=1, barrier=(78.0541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0296167,'amu*angstrom^2'), symmetry=1, barrier=(10.7735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.468516,'amu*angstrom^2'), symmetry=1, barrier=(10.7721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0154947,'amu*angstrom^2'), symmetry=1, barrier=(78.0514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341731,0.0710699,-4.85002e-05,1.83871e-08,-2.91058e-12,36343.3,32.1028], Tmin=(100,'K'), Tmax=(1467.24,'K')), NASAPolynomial(coeffs=[13.2106,0.035987,-1.26341e-05,2.09074e-09,-1.33892e-13,32566.9,-34.9156], Tmin=(1467.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH]C(623)',
    structure = SMILES('[CH2][CH]CC[CH][CH]C'),
    E0 = (572.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,2091.51,2338.86,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36111,0.0497551,-1.77213e-05,1.52484e-09,1.23143e-13,68896.7,28.2554], Tmin=(100,'K'), Tmax=(2624.33,'K')), NASAPolynomial(coeffs=[40.8399,0.00972266,-4.60954e-06,6.75854e-10,-3.24064e-14,43339.4,-199.933], Tmin=(2624.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH][CH]C(555)',
    structure = SMILES('[CH2][CH]C[CH][CH]C'),
    E0 = (596.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,3607.59,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79917,0.0586099,-7.80163e-05,8.14675e-08,-3.30993e-11,71839.7,31.8124], Tmin=(100,'K'), Tmax=(852.592,'K')), NASAPolynomial(coeffs=[-2.4828,0.0533525,-2.41736e-05,4.49758e-09,-3.05589e-13,73491.1,57.1902], Tmin=(852.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2]C([CH2])[CH2](489)',
    structure = SMILES('[CH2]C([CH2])[CH2]'),
    E0 = (463.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00218901,'amu*angstrom^2'), symmetry=1, barrier=(6.67823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226841,'amu*angstrom^2'), symmetry=1, barrier=(69.1489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290286,'amu*angstrom^2'), symmetry=1, barrier=(6.67424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89364,0.0397791,-3.48191e-05,1.81739e-08,-3.6853e-12,55795.5,20.1671], Tmin=(100,'K'), Tmax=(1417.63,'K')), NASAPolynomial(coeffs=[8.35096,0.0156406,-3.01573e-06,2.72841e-10,-9.09425e-15,54559.4,-11.1417], Tmin=(1417.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH][CH]CC([CH2])[CH2](639)',
    structure = SMILES('[CH][CH]CC([CH2])[CH2]'),
    E0 = (855.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,199.748,1040.96,1367.2,1898.02,2038.31],'cm^-1')),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52842,0.0577316,-4.2467e-05,7.26258e-09,7.53295e-12,102951,29.3345], Tmin=(100,'K'), Tmax=(654.247,'K')), NASAPolynomial(coeffs=[7.16544,0.0340336,-1.28177e-05,2.20259e-09,-1.44655e-13,101983,2.76985], Tmin=(654.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(855.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C([CH2])C[C][CH]C(4125)',
    structure = SMILES('[CH2]C([CH2])C[C][CH]C'),
    E0 = (831.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.58956,0.076596,-7.43173e-05,4.60026e-08,-1.22268e-11,100126,33.8149], Tmin=(100,'K'), Tmax=(896.037,'K')), NASAPolynomial(coeffs=[8.34841,0.0419593,-1.63334e-05,2.86117e-09,-1.89968e-13,98735.7,-2.76531], Tmin=(896.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])C[CH][C]C(4126)',
    structure = SMILES('[CH2]C([CH2])C[CH][C]C'),
    E0 = (831.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621327,0.0771403,-7.97703e-05,5.47698e-08,-1.62584e-11,100125,33.5287], Tmin=(100,'K'), Tmax=(841.202,'K')), NASAPolynomial(coeffs=[7.26316,0.0438986,-1.75365e-05,3.1038e-09,-2.0678e-13,99066.5,2.98313], Tmin=(841.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C[CH][CH]C(4127)',
    structure = SMILES('[CH]C([CH2])C[CH][CH]C'),
    E0 = (820.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,202.301,1102.02,1132.92,1564.72,2182.71,2506.67],'cm^-1')),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845363,0.0752019,-8.54319e-05,6.98602e-08,-2.439e-11,98835,35.4677], Tmin=(100,'K'), Tmax=(839.521,'K')), NASAPolynomial(coeffs=[3.93026,0.0499017,-2.12846e-05,3.87819e-09,-2.61858e-13,98690.6,23.3496], Tmin=(839.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(820.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH][CH]CC1CC1(4128)',
    structure = SMILES('C[CH][CH]CC1CC1'),
    E0 = (326.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27448,0.0518123,-6.91609e-06,-1.40776e-08,5.72021e-12,39411.4,30.7059], Tmin=(100,'K'), Tmax=(1200.17,'K')), NASAPolynomial(coeffs=[8.29291,0.04454,-1.79732e-05,3.25519e-09,-2.21376e-13,36565.8,-9.27102], Tmin=(1200.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1CC([CH]C)C1(4129)',
    structure = SMILES('[CH2]C1CC([CH]C)C1'),
    E0 = (324.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40053,0.0394305,5.01504e-05,-8.63717e-08,3.40844e-11,39086.5,28.3206], Tmin=(100,'K'), Tmax=(979.758,'K')), NASAPolynomial(coeffs=[12.5027,0.0376626,-1.38303e-05,2.54002e-09,-1.81318e-13,34820.4,-35.6829], Tmin=(979.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(C)C[CH][CH]C(3934)',
    structure = SMILES('C=C(C)C[CH][CH]C'),
    E0 = (285.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,224.17,454.882,3495.54],'cm^-1')),
        HinderedRotor(inertia=(0.0315361,'amu*angstrom^2'), symmetry=1, barrier=(1.0536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315361,'amu*angstrom^2'), symmetry=1, barrier=(1.0536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315361,'amu*angstrom^2'), symmetry=1, barrier=(1.0536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315361,'amu*angstrom^2'), symmetry=1, barrier=(1.0536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315361,'amu*angstrom^2'), symmetry=1, barrier=(1.0536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31455,0.0560773,-2.52721e-05,4.8448e-09,-3.43424e-13,34409.9,26.5539], Tmin=(100,'K'), Tmax=(2947.02,'K')), NASAPolynomial(coeffs=[41.0517,0.00850427,-3.60546e-06,5.19722e-10,-2.54089e-14,9404.7,-205.885], Tmin=(2947.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(=C)CC[CH]C(4130)',
    structure = SMILES('[CH2]C(=C)CC[CH]C'),
    E0 = (242.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579527,0.0689518,-4.22274e-05,1.32161e-08,-1.70963e-12,29331,29.935], Tmin=(100,'K'), Tmax=(1730.82,'K')), NASAPolynomial(coeffs=[14.4323,0.0369372,-1.44821e-05,2.5293e-09,-1.66019e-13,24535.7,-44.4962], Tmin=(1730.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C=C[CH]C(3935)',
    structure = SMILES('[CH2]C(C)[CH]C=CC'),
    E0 = (237.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,816.323],'cm^-1')),
        HinderedRotor(inertia=(0.862611,'amu*angstrom^2'), symmetry=1, barrier=(19.8331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768887,'amu*angstrom^2'), symmetry=1, barrier=(1.76782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.14136,'amu*angstrom^2'), symmetry=1, barrier=(95.2181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0769576,'amu*angstrom^2'), symmetry=1, barrier=(1.76941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00838642,'amu*angstrom^2'), symmetry=1, barrier=(95.2199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738257,0.0617201,-1.90456e-05,-1.21925e-08,7.42917e-12,28636,28.4776], Tmin=(100,'K'), Tmax=(1058.63,'K')), NASAPolynomial(coeffs=[11.722,0.0396306,-1.5252e-05,2.74009e-09,-1.87796e-13,25222.7,-30.2758], Tmin=(1058.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C[CH][CH]C(4131)',
    structure = SMILES('[CH2]C(=C)C[CH][CH]C'),
    E0 = (437.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,353.954,504.466,3122.48],'cm^-1')),
        HinderedRotor(inertia=(0.0262101,'amu*angstrom^2'), symmetry=1, barrier=(2.54416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0262101,'amu*angstrom^2'), symmetry=1, barrier=(2.54416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0262101,'amu*angstrom^2'), symmetry=1, barrier=(2.54416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0262101,'amu*angstrom^2'), symmetry=1, barrier=(2.54416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0262101,'amu*angstrom^2'), symmetry=1, barrier=(2.54416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73751,0.0569471,-1.10477e-05,-4.75668e-08,4.29254e-11,52665.5,28.4448], Tmin=(100,'K'), Tmax=(511.763,'K')), NASAPolynomial(coeffs=[3.67747,0.0513688,-2.27909e-05,4.32783e-09,-3.03421e-13,52341.4,19.1589], Tmin=(511.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C=C[CH]C(4132)',
    structure = SMILES('[CH2]C([CH2])[CH]C=CC'),
    E0 = (442.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,370.142,1828.5],'cm^-1')),
        HinderedRotor(inertia=(0.0863706,'amu*angstrom^2'), symmetry=1, barrier=(8.39721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00353927,'amu*angstrom^2'), symmetry=1, barrier=(8.39721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00123045,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742144,'amu*angstrom^2'), symmetry=1, barrier=(72.1542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.74217,'amu*angstrom^2'), symmetry=1, barrier=(72.154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.77475,0.063227,-3.11219e-05,5.55221e-11,3.78964e-12,53298.5,29.4327], Tmin=(100,'K'), Tmax=(1020.66,'K')), NASAPolynomial(coeffs=[11.2087,0.0367539,-1.34047e-05,2.32307e-09,-1.55699e-13,50417.6,-24.7971], Tmin=(1020.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C[CH]C=C(4133)',
    structure = SMILES('[CH2]C=CCC([CH2])[CH2]'),
    E0 = (452.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1851.07],'cm^-1')),
        HinderedRotor(inertia=(0.042457,'amu*angstrom^2'), symmetry=1, barrier=(7.52223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186723,'amu*angstrom^2'), symmetry=1, barrier=(3.30858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413131,'amu*angstrom^2'), symmetry=1, barrier=(73.1992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0424598,'amu*angstrom^2'), symmetry=1, barrier=(7.5221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413002,'amu*angstrom^2'), symmetry=1, barrier=(73.1992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621775,0.0654198,-2.96318e-05,-7.28979e-09,7.89136e-12,54554.5,31.3689], Tmin=(100,'K'), Tmax=(964.223,'K')), NASAPolynomial(coeffs=[13.1013,0.0337872,-1.17494e-05,2.00594e-09,-1.34639e-13,51211.8,-33.237], Tmin=(964.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=CC(553)',
    structure = SMILES('[CH2][CH]CC=CC'),
    E0 = (323.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1782.86],'cm^-1')),
        HinderedRotor(inertia=(0.133732,'amu*angstrom^2'), symmetry=1, barrier=(3.07475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133557,'amu*angstrom^2'), symmetry=1, barrier=(3.07073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134216,'amu*angstrom^2'), symmetry=1, barrier=(3.0859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13443,'amu*angstrom^2'), symmetry=1, barrier=(3.09082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80558,0.0490136,-2.40298e-05,5.39389e-09,-4.72015e-13,38969.5,26.5298], Tmin=(100,'K'), Tmax=(2579.8,'K')), NASAPolynomial(coeffs=[18.9745,0.022393,-8.55148e-06,1.39402e-09,-8.43992e-14,30111.1,-72.5714], Tmin=(2579.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][C]([CH2])C[CH][CH]C(4134)',
    structure = SMILES('[CH2][C]([CH2])C[CH][CH]C'),
    E0 = (763.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1881.4,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953113,0.0818466,-0.00012743,1.27439e-07,-4.77692e-11,91881.7,36.2985], Tmin=(100,'K'), Tmax=(892.802,'K')), NASAPolynomial(coeffs=[-2.06846,0.0582179,-2.52885e-05,4.54158e-09,-2.99396e-13,93902.5,58.8287], Tmin=(892.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH][CH][CH]C(4135)',
    structure = SMILES('[CH2]C([CH2])[CH][CH][CH]C'),
    E0 = (772.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,470.946,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06545,0.0731073,-9.15492e-05,8.306e-08,-3.03861e-11,92980.3,37.3174], Tmin=(100,'K'), Tmax=(874.672,'K')), NASAPolynomial(coeffs=[1.06209,0.0529653,-2.24387e-05,4.03682e-09,-2.68864e-13,93752,41.741], Tmin=(874.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898412,0.0752073,-9.10451e-05,7.75073e-08,-2.6973e-11,94275.3,37.1477], Tmin=(100,'K'), Tmax=(883.59,'K')), NASAPolynomial(coeffs=[2.98448,0.0496822,-2.0413e-05,3.61754e-09,-2.38926e-13,94534.4,30.8942], Tmin=(883.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C)C[CH][CH]C(3937)',
    structure = SMILES('[CH2][C](C)C[CH][CH]C'),
    E0 = (558.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,302.098,3682.6,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00116595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116595,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.99347,0.0480048,-1.65089e-05,1.35934e-09,1.03615e-13,67069.7,23.9782], Tmin=(100,'K'), Tmax=(2822.88,'K')), NASAPolynomial(coeffs=[56.4168,-0.00620893,1.63421e-06,-4.06834e-10,3.6978e-14,29476.3,-297.501], Tmin=(2822.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH]C[CH]C(4137)',
    structure = SMILES('[CH2]C([CH2])[CH]C[CH]C'),
    E0 = (577.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,336.471,1289.71,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00140924,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140924,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140924,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140924,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140924,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140924,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54607,0.0615804,-1.36546e-05,-5.71179e-08,5.39183e-11,69575.9,33.1602], Tmin=(100,'K'), Tmax=(521.355,'K')), NASAPolynomial(coeffs=[4.93354,0.0497645,-2.04386e-05,3.70277e-09,-2.5147e-13,69030,17.1766], Tmin=(521.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])[CH2](608)',
    structure = SMILES('[CH2]C[CH]CC([CH2])[CH2]'),
    E0 = (588.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,281.092,942.103,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00200192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00200192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00200192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00200192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00200192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00200192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27897,0.0665019,-3.16281e-05,-2.44679e-08,3.06994e-11,70875.6,33.1399], Tmin=(100,'K'), Tmax=(566.516,'K')), NASAPolynomial(coeffs=[5.74837,0.0484592,-1.96383e-05,3.53139e-09,-2.3878e-13,70152.3,12.2034], Tmin=(566.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)[CH][CH][CH]C(3938)',
    structure = SMILES('[CH2]C(C)[CH][CH][CH]C'),
    E0 = (567.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1230.37,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0133954,'amu*angstrom^2'), symmetry=1, barrier=(14.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133954,'amu*angstrom^2'), symmetry=1, barrier=(14.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133954,'amu*angstrom^2'), symmetry=1, barrier=(14.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133954,'amu*angstrom^2'), symmetry=1, barrier=(14.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133954,'amu*angstrom^2'), symmetry=1, barrier=(14.0346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133954,'amu*angstrom^2'), symmetry=1, barrier=(14.0346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79667,0.0529393,-2.19233e-05,3.64015e-09,-1.99574e-13,68231.9,29.8374], Tmin=(100,'K'), Tmax=(2890.31,'K')), NASAPolynomial(coeffs=[44.5732,0.00487783,-2.04291e-06,2.22285e-10,-4.93993e-15,40008.1,-223.099], Tmin=(2890.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])CC[CH]C(4138)',
    structure = SMILES('[CH2][C]([CH2])CC[CH]C'),
    E0 = (568.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.691666,0.081458,-0.000103783,9.16878e-08,-3.26421e-11,68508.9,34.6736], Tmin=(100,'K'), Tmax=(878.136,'K')), NASAPolynomial(coeffs=[2.2729,0.0541354,-2.27433e-05,4.07219e-09,-2.70352e-13,69007,31.6675], Tmin=(878.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH][CH]CC(606)',
    structure = SMILES('[CH2]C([CH2])[CH][CH]CC'),
    E0 = (577.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,350.462,3118.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0231768,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231768,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231768,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231768,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231768,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231768,'amu*angstrom^2'), symmetry=1, barrier=(3.66408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10946,0.0522036,4.38182e-05,-2.09676e-07,1.92268e-10,69553.1,31.1501], Tmin=(100,'K'), Tmax=(420.994,'K')), NASAPolynomial(coeffs=[4.32374,0.0507899,-2.10676e-05,3.80126e-09,-2.55774e-13,69192.7,20.3177], Tmin=(420.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])[CH2](4139)',
    structure = SMILES('[CH2][CH]CCC([CH2])[CH2]'),
    E0 = (588.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,376.844,1031.52,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0372877,'amu*angstrom^2'), symmetry=1, barrier=(1.13787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372877,'amu*angstrom^2'), symmetry=1, barrier=(1.13787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372877,'amu*angstrom^2'), symmetry=1, barrier=(1.13787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372877,'amu*angstrom^2'), symmetry=1, barrier=(1.13787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372877,'amu*angstrom^2'), symmetry=1, barrier=(1.13787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372877,'amu*angstrom^2'), symmetry=1, barrier=(1.13787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723893,0.0737262,-6.32404e-05,3.58314e-08,-9.05529e-12,70898.9,35.2153], Tmin=(100,'K'), Tmax=(922.801,'K')), NASAPolynomial(coeffs=[7.13554,0.0459343,-1.80655e-05,3.19562e-09,-2.13868e-13,69715.5,4.79786], Tmin=(922.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH2])C[CH]CC(605)',
    structure = SMILES('[CH2][C]([CH2])C[CH]CC'),
    E0 = (568.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300.704,2672.4,3846.6],'cm^-1')),
        HinderedRotor(inertia=(0.002775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.002775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.002775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.002775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.002775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.002775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735946,0.0818459,-0.000108642,9.9604e-08,-3.62669e-11,68507.5,34.3431], Tmin=(100,'K'), Tmax=(876.766,'K')), NASAPolynomial(coeffs=[1.17352,0.056099,-2.39606e-05,4.31818e-09,-2.87443e-13,69343.6,37.4955], Tmin=(876.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])C(3939)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])C'),
    E0 = (577.894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,354.469,3058.14,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0247994,'amu*angstrom^2'), symmetry=1, barrier=(3.74496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247994,'amu*angstrom^2'), symmetry=1, barrier=(3.74496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247994,'amu*angstrom^2'), symmetry=1, barrier=(3.74496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247994,'amu*angstrom^2'), symmetry=1, barrier=(3.74496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247994,'amu*angstrom^2'), symmetry=1, barrier=(3.74496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247994,'amu*angstrom^2'), symmetry=1, barrier=(3.74496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954522,0.0726474,-7.55607e-05,6.14412e-08,-2.20725e-11,69608.8,35.8587], Tmin=(100,'K'), Tmax=(827.656,'K')), NASAPolynomial(coeffs=[2.73208,0.053852,-2.30026e-05,4.20939e-09,-2.85491e-13,69664.1,29.7308], Tmin=(827.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (577.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (577.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (735.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (743.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1034.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1041.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1047.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1043.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1043.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1032.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (585.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (586.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (600.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (641.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (641.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (660.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (661.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (672.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (705.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (627.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (643.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (934.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (974.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (984.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (994.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (719.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (714.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (730.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (695.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (693.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (724.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (729.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (695.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (630.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['allyl(82)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2]C([CH2])CC=CC(604)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2][CH]CC[CH][CH]C(623)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH2])C([CH2])[CH]C(570)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH2][CH]C[CH][CH]C(555)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH]C(3874)', '[CH2]C([CH2])[CH2](489)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.02491e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH][CH]CC([CH2])[CH2](639)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2]C([CH2])C[C][CH]C(4125)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2]C([CH2])C[CH][C]C(4126)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([CH2])C[CH][CH]C(4127)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['C[CH][CH]CC1CC1(4128)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2]C1CC([CH]C)C1(4129)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['C=C(C)C[CH][CH]C(3934)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2]C(=C)CC[CH]C(4130)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2]C(C)C=C[CH]C(3935)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C(=C)C[CH][CH]C(4131)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]C([CH2])C=C[CH]C(4132)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C([CH2])C[CH]C=C(4133)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(T)(20)', '[CH2][CH]CC=CC(553)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['allyl(82)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][CH2](497)', 'm1_allyl(186)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C]([CH2])C[CH][CH]C(4134)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C([CH2])[CH][CH][CH]C(4135)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2][C](C)C[CH][CH]C(3937)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])[CH]C[CH]C(4137)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C[CH]CC([CH2])[CH2](608)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2]C(C)[CH][CH][CH]C(3938)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2][C]([CH2])CC[CH]C(4138)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH2])[CH][CH]CC(606)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]CCC([CH2])[CH2](4139)'],
    products = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2][C]([CH2])C[CH]CC(605)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(10873.3,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1110',
    isomers = [
        '[CH2]C([CH2])C[CH][CH]C(607)',
    ],
    reactants = [
        ('allyl(82)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1110',
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


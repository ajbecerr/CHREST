species(
    label = '[CH2][CH]CC(C)[CH][CH2](1276)',
    structure = SMILES('[CH2][CH]CC(C)[CH][CH2]'),
    E0 = (581.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,363.316,3103.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28585,0.0653206,-4.33609e-05,1.8603e-08,-4.13437e-12,70030,34.6848], Tmin=(100,'K'), Tmax=(919.114,'K')), NASAPolynomial(coeffs=[4.00684,0.0534791,-2.40359e-05,4.58607e-09,-3.2183e-13,69529.8,21.7871], Tmin=(919.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH]CC(C)C=C(1273)',
    structure = SMILES('[CH2][CH]CC(C)C=C'),
    E0 = (303.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.68982,0.0656371,-3.76021e-05,1.07496e-08,-1.25464e-12,36670.6,32.9921], Tmin=(100,'K'), Tmax=(1912.45,'K')), NASAPolynomial(coeffs=[15.5232,0.034612,-1.3268e-05,2.26686e-09,-1.45749e-13,30997,-48.1884], Tmin=(1912.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]CC[CH][CH2](483)',
    structure = SMILES('[CH2][CH]CC[CH][CH2]'),
    E0 = (607.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,406.91,2706.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.63,'J/mol'), sigma=(6.0508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.24 K, Pc=32.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96636,0.0441879,-1.8776e-05,3.20919e-09,-1.84564e-13,73082,26.0155], Tmin=(100,'K'), Tmax=(2870.12,'K')), NASAPolynomial(coeffs=[37.1311,0.00432169,-1.9902e-06,2.50805e-10,-8.8047e-15,50279.3,-180.391], Tmin=(2870.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])[CH]C(586)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH]C'),
    E0 = (581.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,353.512,3195.31,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5058,0.0619309,-2.04238e-05,-3.53119e-08,3.53303e-11,70000.3,33.9722], Tmin=(100,'K'), Tmax=(536.631,'K')), NASAPolynomial(coeffs=[4.26044,0.0520363,-2.25022e-05,4.21175e-09,-2.92542e-13,69551.5,20.9701], Tmin=(536.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])C(4096)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])C'),
    E0 = (577.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5058,0.0619309,-2.04238e-05,-3.53119e-08,3.53303e-11,69597.8,33.9722], Tmin=(100,'K'), Tmax=(536.631,'K')), NASAPolynomial(coeffs=[4.26044,0.0520363,-2.25022e-05,4.21175e-09,-2.92542e-13,69148.9,20.9701], Tmin=(536.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])[CH2](1265)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])[CH2]'),
    E0 = (586.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,341.472,3095.34],'cm^-1')),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3395.46,'J/mol'), sigma=(6.43099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.36 K, Pc=28.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596231,0.0714284,-5.11594e-05,2.10748e-08,-3.72485e-12,70636.9,35.2063], Tmin=(100,'K'), Tmax=(1296.67,'K')), NASAPolynomial(coeffs=[10.6242,0.0404935,-1.53733e-05,2.67563e-09,-1.77439e-13,68036.3,-15.7783], Tmin=(1296.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH]C(124)',
    structure = SMILES('[CH2][CH]C[CH]C'),
    E0 = (426.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1955.55,1958.71],'cm^-1')),
        HinderedRotor(inertia=(0.00224622,'amu*angstrom^2'), symmetry=1, barrier=(6.10841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.024144,'amu*angstrom^2'), symmetry=1, barrier=(65.4821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120719,'amu*angstrom^2'), symmetry=1, barrier=(6.15289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120526,'amu*angstrom^2'), symmetry=1, barrier=(6.14858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95832,0.0347404,-1.13952e-05,5.01961e-10,2.04138e-13,51268.6,22.8061], Tmin=(100,'K'), Tmax=(2290.84,'K')), NASAPolynomial(coeffs=[16.7487,0.0200393,-7.90976e-06,1.27465e-09,-7.52019e-14,42489.5,-60.5265], Tmin=(2290.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][C]CC(C)[CH][CH2](4097)',
    structure = SMILES('[CH2][C]CC(C)[CH][CH2]'),
    E0 = (835.258,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938742,0.0696321,-4.98365e-05,1.93416e-08,-3.24935e-12,100566,33.0753], Tmin=(100,'K'), Tmax=(1319.94,'K')), NASAPolynomial(coeffs=[10.0782,0.0419357,-1.83622e-05,3.44493e-09,-2.38505e-13,98153.5,-13.5544], Tmin=(1319.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]C(C)C[CH][CH2](4098)',
    structure = SMILES('[CH2][C]C(C)C[CH][CH2]'),
    E0 = (835.162,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.957224,0.0707496,-5.61645e-05,2.69475e-08,-5.82617e-12,100553,33.2032], Tmin=(100,'K'), Tmax=(1042.78,'K')), NASAPolynomial(coeffs=[7.54038,0.0454976,-1.98408e-05,3.72544e-09,-2.58895e-13,99180.2,1.16741], Tmin=(1042.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(835.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC(C)[CH][CH2](4099)',
    structure = SMILES('[CH][CH]CC(C)[CH][CH2]'),
    E0 = (824.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,219.891,1243.09,1313.89,1415.42,1703.6,1765.87],'cm^-1')),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11937,0.0685479,-5.55915e-05,3.01639e-08,-7.86861e-12,99258.6,34.4997], Tmin=(100,'K'), Tmax=(855.322,'K')), NASAPolynomial(coeffs=[5.17,0.0496042,-2.23685e-05,4.26817e-09,-2.99385e-13,98565.7,15.5909], Tmin=(855.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C(C)C[CH][CH2](4100)',
    structure = SMILES('[CH][CH]C(C)C[CH][CH2]'),
    E0 = (824.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,219.891,1243.09,1313.89,1415.42,1703.6,1765.87],'cm^-1')),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102584,'amu*angstrom^2'), symmetry=1, barrier=(3.09541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11937,0.0685479,-5.55915e-05,3.01639e-08,-7.86861e-12,99258.6,34.4997], Tmin=(100,'K'), Tmax=(855.322,'K')), NASAPolynomial(coeffs=[5.17,0.0496042,-2.23685e-05,4.26817e-09,-2.99385e-13,98565.7,15.5909], Tmin=(855.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C=C(C)CC[CH2](4102)',
    structure = SMILES('[CH2]CC[C](C)C=C'),
    E0 = (241.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778085,0.0602103,-1.46318e-05,-1.62658e-08,8.59457e-12,29157.9,28.1738], Tmin=(100,'K'), Tmax=(1068.99,'K')), NASAPolynomial(coeffs=[11.8681,0.039827,-1.56567e-05,2.84954e-09,-1.96793e-13,25580.5,-31.7114], Tmin=(1068.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC(C)C[CH2](4103)',
    structure = SMILES('[CH2]C=CC(C)C[CH2]'),
    E0 = (248.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.577771,0.0629928,-1.50837e-05,-2.06048e-08,1.10905e-11,30041.4,30.4804], Tmin=(100,'K'), Tmax=(1040.01,'K')), NASAPolynomial(coeffs=[13.9787,0.0370218,-1.45065e-05,2.66636e-09,-1.86422e-13,25871.1,-41.3456], Tmin=(1040.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C[CH][CH2](653)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH2]'),
    E0 = (801.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2055.29,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82635,0.0594615,-8.88659e-05,9.41844e-08,-3.76943e-11,96522.8,33.5013], Tmin=(100,'K'), Tmax=(862.635,'K')), NASAPolynomial(coeffs=[-2.63202,0.0509582,-2.33462e-05,4.3408e-09,-2.93791e-13,98377.6,60.6438], Tmin=(862.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(801.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C](C)[CH][CH2](4106)',
    structure = SMILES('[CH2][CH]C[C](C)[CH][CH2]'),
    E0 = (766.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,3732.07],'cm^-1')),
        HinderedRotor(inertia=(0.00533923,'amu*angstrom^2'), symmetry=1, barrier=(6.51025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00533923,'amu*angstrom^2'), symmetry=1, barrier=(6.51025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00533923,'amu*angstrom^2'), symmetry=1, barrier=(6.51025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00533923,'amu*angstrom^2'), symmetry=1, barrier=(6.51025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00533923,'amu*angstrom^2'), symmetry=1, barrier=(6.51025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00533923,'amu*angstrom^2'), symmetry=1, barrier=(6.51025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.69282,0.0489295,-1.99249e-05,2.71222e-09,-6.04585e-14,92204.6,26.8919], Tmin=(100,'K'), Tmax=(2832.23,'K')), NASAPolynomial(coeffs=[54.5607,-0.00683832,1.09789e-06,-2.3243e-10,2.25888e-14,56943.9,-282.854], Tmin=(2832.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(766.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C(C)[CH][CH2](4107)',
    structure = SMILES('[CH2][CH][CH]C(C)[CH][CH2]'),
    E0 = (776.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3207.06,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0243687,'amu*angstrom^2'), symmetry=1, barrier=(6.23651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0243687,'amu*angstrom^2'), symmetry=1, barrier=(6.23651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0243687,'amu*angstrom^2'), symmetry=1, barrier=(6.23651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0243687,'amu*angstrom^2'), symmetry=1, barrier=(6.23651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0243687,'amu*angstrom^2'), symmetry=1, barrier=(6.23651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0243687,'amu*angstrom^2'), symmetry=1, barrier=(6.23651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49734,0.0538544,-2.53221e-05,4.98281e-09,-3.61742e-13,93366.7,32.7458], Tmin=(100,'K'), Tmax=(2976.69,'K')), NASAPolynomial(coeffs=[42.8411,0.00416776,-2.56516e-06,3.96936e-10,-1.94948e-14,67343.2,-209.265], Tmin=(2976.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(776.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[C](C)C[CH2](4109)',
    structure = SMILES('[CH2][CH]C[C](C)C[CH2]'),
    E0 = (572.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.814567,0.0786134,-9.68146e-05,8.84874e-08,-3.3797e-11,68946.6,34.9334], Tmin=(100,'K'), Tmax=(819.511,'K')), NASAPolynomial(coeffs=[1.53158,0.0575737,-2.61996e-05,4.92565e-09,-3.38466e-13,69418.1,35.2105], Tmin=(819.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C[CH2](4110)',
    structure = SMILES('[CH2][CH]C(C)[CH]C[CH2]'),
    E0 = (581.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,2049.5,3553.19],'cm^-1')),
        HinderedRotor(inertia=(0.0933973,'amu*angstrom^2'), symmetry=1, barrier=(2.14739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933973,'amu*angstrom^2'), symmetry=1, barrier=(2.14739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933973,'amu*angstrom^2'), symmetry=1, barrier=(2.14739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933973,'amu*angstrom^2'), symmetry=1, barrier=(2.14739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933973,'amu*angstrom^2'), symmetry=1, barrier=(2.14739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933973,'amu*angstrom^2'), symmetry=1, barrier=(2.14739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43023,0.0628182,-3.39647e-05,8.52602e-09,-8.5982e-13,70033.5,33.926], Tmin=(100,'K'), Tmax=(2091.1,'K')), NASAPolynomial(coeffs=[14.1144,0.0385551,-1.65602e-05,2.97726e-09,-1.96441e-13,64728.7,-36.6245], Tmin=(2091.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C(C)C[CH2](4111)',
    structure = SMILES('[CH2][CH][CH]C(C)C[CH2]'),
    E0 = (581.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,363.316,3103.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0269495,'amu*angstrom^2'), symmetry=1, barrier=(4.00742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28585,0.0653206,-4.33609e-05,1.8603e-08,-4.13437e-12,70030,34.6848], Tmin=(100,'K'), Tmax=(919.114,'K')), NASAPolynomial(coeffs=[4.00684,0.0534791,-2.40359e-05,4.58607e-09,-3.2183e-13,69529.8,21.7871], Tmin=(919.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])C[CH2](588)',
    structure = SMILES('[CH2][CH]CC([CH2])C[CH2]'),
    E0 = (592.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,181.491,1620.2,3606.77],'cm^-1')),
        HinderedRotor(inertia=(0.0566294,'amu*angstrom^2'), symmetry=1, barrier=(2.29896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0566294,'amu*angstrom^2'), symmetry=1, barrier=(2.29896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0566294,'amu*angstrom^2'), symmetry=1, barrier=(2.29896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0566294,'amu*angstrom^2'), symmetry=1, barrier=(2.29896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0566294,'amu*angstrom^2'), symmetry=1, barrier=(2.29896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0566294,'amu*angstrom^2'), symmetry=1, barrier=(2.29896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907477,0.0708075,-5.41887e-05,2.63505e-08,-5.88168e-12,71313.7,35.2568], Tmin=(100,'K'), Tmax=(1008.42,'K')), NASAPolynomial(coeffs=[6.73307,0.0476997,-1.98165e-05,3.62704e-09,-2.48266e-13,70138.7,7.10285], Tmin=(1008.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C](C)CC[CH2](4112)',
    structure = SMILES('[CH2][CH][C](C)CC[CH2]'),
    E0 = (572.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.830422,0.0769041,-8.78214e-05,7.71414e-08,-2.97532e-11,68958.9,34.698], Tmin=(100,'K'), Tmax=(778.736,'K')), NASAPolynomial(coeffs=[2.44559,0.0567962,-2.63383e-05,5.02939e-09,-3.50095e-13,69065.5,29.6092], Tmin=(778.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C](C)[CH]C(1274)',
    structure = SMILES('[CH2][CH]C[C](C)[CH]C'),
    E0 = (561.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,417.121,1262.57,4000],'cm^-1')),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24131,0.0523227,-2.11878e-05,2.88029e-09,-5.64948e-14,67543,26.7888], Tmin=(100,'K'), Tmax=(2714.8,'K')), NASAPolynomial(coeffs=[46.2883,0.00396482,-2.7941e-06,4.07837e-10,-1.70804e-14,38617.7,-234.106], Tmin=(2714.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][CH]C(3909)',
    structure = SMILES('[CH2][CH]C(C)[CH][CH]C'),
    E0 = (570.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,549.259,3239.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04461,0.0572562,-2.65999e-05,5.15936e-09,-3.5928e-13,68705.3,32.6477], Tmin=(100,'K'), Tmax=(2641.81,'K')), NASAPolynomial(coeffs=[33.6327,0.0159202,-6.81576e-06,1.09699e-09,-6.28765e-14,49749.9,-154.721], Tmin=(2641.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])CC[CH2](4113)',
    structure = SMILES('[CH2][CH]C([CH2])CC[CH2]'),
    E0 = (592.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,228.891,1037.49,3344.88],'cm^-1')),
        HinderedRotor(inertia=(0.0471359,'amu*angstrom^2'), symmetry=1, barrier=(1.52149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0471359,'amu*angstrom^2'), symmetry=1, barrier=(1.52149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0471359,'amu*angstrom^2'), symmetry=1, barrier=(1.52149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0471359,'amu*angstrom^2'), symmetry=1, barrier=(1.52149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0471359,'amu*angstrom^2'), symmetry=1, barrier=(1.52149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0471359,'amu*angstrom^2'), symmetry=1, barrier=(1.52149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.893968,0.0696099,-4.75125e-05,1.82446e-08,-3.08601e-12,71326.7,35.113], Tmin=(100,'K'), Tmax=(1300.91,'K')), NASAPolynomial(coeffs=[9.06999,0.0444708,-1.85264e-05,3.39049e-09,-2.31487e-13,69199.4,-6.48252], Tmin=(1300.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C(C)[CH]C(1275)',
    structure = SMILES('[CH2][CH][CH]C(C)[CH]C'),
    E0 = (570.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,549.259,3239.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04461,0.0572562,-2.65999e-05,5.15936e-09,-3.5928e-13,68705.3,32.6477], Tmin=(100,'K'), Tmax=(2641.81,'K')), NASAPolynomial(coeffs=[33.6327,0.0159202,-6.81576e-06,1.09699e-09,-6.28765e-14,49749.9,-154.721], Tmin=(2641.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C](C)C[CH]C(3908)',
    structure = SMILES('[CH2][CH][C](C)C[CH]C'),
    E0 = (561.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,417.121,1262.57,4000],'cm^-1')),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012739,'amu*angstrom^2'), symmetry=1, barrier=(2.70041,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24131,0.0523227,-2.11878e-05,2.88029e-09,-5.64948e-14,67543,26.7888], Tmin=(100,'K'), Tmax=(2714.8,'K')), NASAPolynomial(coeffs=[46.2883,0.00396482,-2.7941e-06,4.07837e-10,-1.70804e-14,38617.7,-234.106], Tmin=(2714.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH]C(3910)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH]C'),
    E0 = (581.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,353.512,3195.31,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50564,0.0619335,-2.04379e-05,-3.52822e-08,3.5309e-11,70000.3,33.9727], Tmin=(100,'K'), Tmax=(536.66,'K')), NASAPolynomial(coeffs=[4.26049,0.0520362,-2.25021e-05,4.21173e-09,-2.92541e-13,69551.5,20.9698], Tmin=(536.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (581.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (581.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1026.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (724.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (776.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (743.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1039.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1044.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1047.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1046.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1036.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1036.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (589.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (644.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (644.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (655.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (662.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (643.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (645.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (643.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (934.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (938.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (978.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (987.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (998.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (747.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (718.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (697.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (709.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (697.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (730.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (722.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (668.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (712.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (712.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (635.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['allyl(82)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]CC(C)C=C(1273)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]CC[CH][CH2](483)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(287528,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]CC([CH2])[CH]C(586)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]C[CH]C([CH2])C(4096)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]C(C)C([CH2])[CH2](1265)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C[CH]C(124)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C([CH2])C(111)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]CC(C)[CH][CH2](4097)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][C]C(C)C[CH][CH2](4098)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH]CC(C)[CH][CH2](4099)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH][CH]C(C)C[CH][CH2](4100)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2]C1CC(C)C1[CH2](4101)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2]C=C(C)CC[CH2](4102)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2]C=CC(C)C[CH2](4103)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2][CH]CC(C)=C[CH2](4104)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2][CH]C(C)C=C[CH2](4105)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH][CH2](497)', 'm1_allyl(186)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH3(17)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.139979,'m^3/(mol*s)'), n=2.09962, Ea=(33.817,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['allyl(82)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3(17)', '[CH2][CH][CH]C[CH][CH2](653)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][CH]C[C](C)[CH][CH2](4106)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH][CH]C(C)[CH][CH2](4107)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH]CC([CH2])[CH][CH2](4108)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]C[C](C)C[CH2](4109)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]C(C)[CH]C[CH2](4110)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH]C(C)C[CH2](4111)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]CC([CH2])C[CH2](588)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH][C](C)CC[CH2](4112)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]C[C](C)[CH]C(1274)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.081e+09,'s^-1'), n=0.812344, Ea=(148.932,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_Cs2] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]C(C)[CH][CH]C(3909)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([CH2])CC[CH2](4113)'],
    products = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH][CH]C(C)[CH]C(1275)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH][C](C)C[CH]C(3908)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]CC(C)[CH][CH2](1276)'],
    products = ['[CH2][CH]C([CH2])C[CH]C(3910)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1117',
    isomers = [
        '[CH2][CH]CC(C)[CH][CH2](1276)',
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
    label = 'PDepNetwork #1117',
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


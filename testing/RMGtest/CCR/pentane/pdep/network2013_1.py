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
    label = '[CH2][CH][CH]CC[CH2](655)',
    structure = SMILES('[CH2][CH][CH]CC[CH2]'),
    E0 = (607.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,502.332,2610.7,3566.25],'cm^-1')),
        HinderedRotor(inertia=(0.0124219,'amu*angstrom^2'), symmetry=1, barrier=(2.01271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0124219,'amu*angstrom^2'), symmetry=1, barrier=(2.01271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0124219,'amu*angstrom^2'), symmetry=1, barrier=(2.01271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0124219,'amu*angstrom^2'), symmetry=1, barrier=(2.01271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0124219,'amu*angstrom^2'), symmetry=1, barrier=(2.01271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.63,'J/mol'), sigma=(6.0508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.24 K, Pc=32.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18168,0.0425548,-1.66264e-05,2.15789e-09,-2.55685e-14,73073.4,25.7665], Tmin=(100,'K'), Tmax=(2657.42,'K')), NASAPolynomial(coeffs=[33.3343,0.00857218,-3.88157e-06,5.75411e-10,-2.86112e-14,53021.2,-156.748], Tmin=(2657.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH]C(3964)',
    structure = SMILES('[CH2][CH][CH]CC[CH]C'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36111,0.0497551,-1.77213e-05,1.52485e-09,1.23143e-13,68896.7,28.2554], Tmin=(100,'K'), Tmax=(2624.33,'K')), NASAPolynomial(coeffs=[40.8399,0.00972263,-4.60953e-06,6.75852e-10,-3.24062e-14,43339.4,-199.933], Tmin=(2624.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]CC(622)',
    structure = SMILES('[CH2][CH][CH]C[CH]CC'),
    E0 = (572.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,2351.73,2410.93,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15593,0.0742926,-9.87011e-05,1.0042e-07,-4.03765e-11,69002.7,36.1938], Tmin=(100,'K'), Tmax=(843.06,'K')), NASAPolynomial(coeffs=[-2.56545,0.0635643,-2.91099e-05,5.4533e-09,-3.72438e-13,70638.9,59.4945], Tmin=(843.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])C(3893)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])C'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596231,0.0714284,-5.11594e-05,2.10748e-08,-3.72485e-12,70636.9,35.8994], Tmin=(100,'K'), Tmax=(1296.67,'K')), NASAPolynomial(coeffs=[10.6242,0.0404935,-1.53733e-05,2.67563e-09,-1.77439e-13,68036.3,-15.0851], Tmin=(1296.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[CH]C(556)',
    structure = SMILES('[CH2][CH][CH]C[CH]C'),
    E0 = (596.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,3630.26,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79918,0.0586097,-7.80157e-05,8.14666e-08,-3.30989e-11,71839.7,31.8124], Tmin=(100,'K'), Tmax=(852.594,'K')), NASAPolynomial(coeffs=[-2.48283,0.0533526,-2.41736e-05,4.49758e-09,-3.0559e-13,73491.1,57.1903], Tmin=(852.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH][CH2](5537)',
    structure = SMILES('[CH][CH][CH2]'),
    E0 = (727.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1609,1609.01,1609.02],'cm^-1')),
        HinderedRotor(inertia=(0.0337841,'amu*angstrom^2'), symmetry=1, barrier=(5.78124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00314687,'amu*angstrom^2'), symmetry=1, barrier=(5.78111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42424,0.0149081,-5.18362e-06,2.15834e-10,1.32851e-13,87550.9,15.303], Tmin=(100,'K'), Tmax=(1973.84,'K')), NASAPolynomial(coeffs=[8.10291,0.00913487,-3.61425e-06,6.37549e-10,-4.11102e-14,84981.5,-12.2801], Tmin=(1973.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])C(102)',
    structure = SMILES('[CH2]C([CH2])C'),
    E0 = (258.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.225531,'amu*angstrom^2'), symmetry=1, barrier=(5.39243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211324,'amu*angstrom^2'), symmetry=1, barrier=(5.37273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00963902,'amu*angstrom^2'), symmetry=1, barrier=(80.4753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2657,0.0333191,-4.86449e-06,-1.80964e-08,1.05994e-11,31115.8,18.163], Tmin=(100,'K'), Tmax=(878.748,'K')), NASAPolynomial(coeffs=[8.49826,0.0192676,-5.32045e-06,7.92236e-10,-4.96039e-14,29467.6,-14.2455], Tmin=(878.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2][CH][C]CC([CH2])C(9042)',
    structure = SMILES('[CH2][CH][C]CC([CH2])C'),
    E0 = (831.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76705,0.0737167,-6.52629e-05,3.63482e-08,-8.91782e-12,100139,33.8805], Tmin=(100,'K'), Tmax=(950.059,'K')), NASAPolynomial(coeffs=[7.84809,0.0439039,-1.81935e-05,3.31934e-09,-2.2665e-13,98793.2,0.0812971], Tmin=(950.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC([CH2])C(9043)',
    structure = SMILES('[CH2][C][CH]CC([CH2])C'),
    E0 = (831.663,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794189,0.0742686,-7.05224e-05,4.45657e-08,-1.25911e-11,100138,33.614], Tmin=(100,'K'), Tmax=(829.711,'K')), NASAPolynomial(coeffs=[6.64445,0.0460645,-1.9533e-05,3.59581e-09,-2.46372e-13,99167.2,6.48202], Tmin=(829.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(C)C[CH][CH][CH2](9044)',
    structure = SMILES('[CH]C(C)C[CH][CH][CH2]'),
    E0 = (821.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,195.277,925.785,1038.05,2107.77,2403.71,2755.61],'cm^-1')),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92656,0.073534,-8.10432e-05,6.70608e-08,-2.44693e-11,98851.7,35.1812], Tmin=(100,'K'), Tmax=(795.536,'K')), NASAPolynomial(coeffs=[3.46626,0.051787,-2.31118e-05,4.32889e-09,-2.97939e-13,98731.7,25.295], Tmin=(795.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(821.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC([CH2])C(9045)',
    structure = SMILES('[CH][CH][CH]CC([CH2])C'),
    E0 = (820.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,202.301,1102.02,1132.92,1564.72,2182.71,2506.67],'cm^-1')),
        HinderedRotor(inertia=(0.0623098,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623098,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623098,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623098,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623098,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623098,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845372,0.0752017,-8.54314e-05,6.98595e-08,-2.43897e-11,98835,35.4677], Tmin=(100,'K'), Tmax=(839.525,'K')), NASAPolynomial(coeffs=[3.93024,0.0499017,-2.12847e-05,3.8782e-09,-2.61858e-13,98690.6,23.3497], Tmin=(839.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(820.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CC(C)C1(5009)',
    structure = SMILES('[CH2][CH]C1CC(C)C1'),
    E0 = (324.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41806,0.0384786,5.22954e-05,-8.68703e-08,3.34183e-11,39105.9,28.2638], Tmin=(100,'K'), Tmax=(1000.35,'K')), NASAPolynomial(coeffs=[12.6123,0.0385677,-1.50899e-05,2.85622e-09,-2.06216e-13,34622.2,-36.9624], Tmin=(1000.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C)CCC=C(5027)',
    structure = SMILES('[CH2][C](C)CCC=C'),
    E0 = (296.948,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07964,0.0681199,-4.86681e-05,2.21912e-08,-4.82384e-12,35816.5,29.942], Tmin=(100,'K'), Tmax=(1002.29,'K')), NASAPolynomial(coeffs=[5.67339,0.0497872,-2.12323e-05,3.94274e-09,-2.7222e-13,34895.6,7.76919], Tmin=(1002.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)C(5017)',
    structure = SMILES('[CH2]C=C[CH]C(C)C'),
    E0 = (183.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.758855,0.0570875,3.5624e-06,-4.00003e-08,1.78415e-11,22194.5,26.2091], Tmin=(100,'K'), Tmax=(1019.15,'K')), NASAPolynomial(coeffs=[13.9035,0.0374846,-1.46652e-05,2.71945e-09,-1.92171e-13,17854,-45.6054], Tmin=(1019.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C]([CH2])C(5011)',
    structure = SMILES('[CH2]C=CC[C]([CH2])C'),
    E0 = (432.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.102837,'amu*angstrom^2'), symmetry=1, barrier=(2.36441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995862,'amu*angstrom^2'), symmetry=1, barrier=(2.28968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0984102,'amu*angstrom^2'), symmetry=1, barrier=(2.26264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982723,'amu*angstrom^2'), symmetry=1, barrier=(2.25947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00727311,'amu*angstrom^2'), symmetry=1, barrier=(2.31243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08981,0.0663958,-4.71361e-05,1.99482e-08,-3.82871e-12,52163.1,29.7387], Tmin=(100,'K'), Tmax=(1149.28,'K')), NASAPolynomial(coeffs=[7.35999,0.0445731,-1.86543e-05,3.42689e-09,-2.3492e-13,50721.8,-1.38374], Tmin=(1149.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])C(5012)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])C'),
    E0 = (388.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1451.4,2646.44],'cm^-1')),
        HinderedRotor(inertia=(0.828594,'amu*angstrom^2'), symmetry=1, barrier=(19.051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0371505,'amu*angstrom^2'), symmetry=1, barrier=(19.0516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.48473,'amu*angstrom^2'), symmetry=1, barrier=(80.1207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828569,'amu*angstrom^2'), symmetry=1, barrier=(19.0504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0535999,'amu*angstrom^2'), symmetry=1, barrier=(80.1131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787236,0.0586763,-8.72898e-06,-2.75874e-08,1.41857e-11,46857.4,28.5804], Tmin=(100,'K'), Tmax=(994.507,'K')), NASAPolynomial(coeffs=[13.5045,0.0344226,-1.27148e-05,2.27869e-09,-1.58143e-13,42997.8,-39.3902], Tmin=(994.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C[C]([CH2])C(9046)',
    structure = SMILES('[CH2][CH][CH]C[C]([CH2])C'),
    E0 = (763.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1314.93,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998698,0.0806099,-0.000124564,1.26537e-07,-4.85632e-11,91899.9,36.8318], Tmin=(100,'K'), Tmax=(874.732,'K')), NASAPolynomial(coeffs=[-2.3313,0.0597485,-2.69057e-05,4.94172e-09,-3.31223e-13,93863.1,60.3435], Tmin=(874.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C([CH2])C(9047)',
    structure = SMILES('[CH2][CH][CH][CH]C([CH2])C'),
    E0 = (772.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1692.17,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12398,0.0717174,-8.81604e-05,8.1535e-08,-3.09609e-11,92998,37.8045], Tmin=(100,'K'), Tmax=(845.322,'K')), NASAPolynomial(coeffs=[0.71883,0.0546362,-2.41384e-05,4.4567e-09,-3.02345e-13,93745.3,43.7059], Tmin=(845.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[C](C)C(9048)',
    structure = SMILES('[CH2][CH][CH]C[C](C)C'),
    E0 = (558.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0181,0.0785203,-0.000110911,1.13133e-07,-4.49377e-11,67234.9,34.2858], Tmin=(100,'K'), Tmax=(848.085,'K')), NASAPolynomial(coeffs=[-2.49183,0.0637543,-2.93973e-05,5.50982e-09,-3.75781e-13,68956.7,57.2815], Tmin=(848.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])C(4096)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])C'),
    E0 = (577.978,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5058,0.0619309,-2.04238e-05,-3.53119e-08,3.53303e-11,69597.8,33.9722], Tmin=(100,'K'), Tmax=(536.631,'K')), NASAPolynomial(coeffs=[4.26044,0.0520363,-2.25022e-05,4.21175e-09,-2.92542e-13,69148.9,20.9701], Tmin=(536.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)C(9049)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)C'),
    E0 = (567.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64245,0.0532781,-2.17205e-05,3.10983e-09,-7.84773e-14,68262,29.7664], Tmin=(100,'K'), Tmax=(2563.99,'K')), NASAPolynomial(coeffs=[34.5013,0.0156879,-6.81493e-06,1.07658e-09,-5.98657e-14,47943.8,-161.695], Tmin=(2563.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C]([CH2])C(5005)',
    structure = SMILES('[CH2][CH]CC[C]([CH2])C'),
    E0 = (568.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,473.669,2851.12,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0326955,'amu*angstrom^2'), symmetry=1, barrier=(4.53797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0326955,'amu*angstrom^2'), symmetry=1, barrier=(4.53797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0326955,'amu*angstrom^2'), symmetry=1, barrier=(4.53797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0326955,'amu*angstrom^2'), symmetry=1, barrier=(4.53797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0326955,'amu*angstrom^2'), symmetry=1, barrier=(4.53797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0326955,'amu*angstrom^2'), symmetry=1, barrier=(4.53797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752306,0.0800394,-0.000100271,8.99626e-08,-3.31081e-11,68526.5,35.1533], Tmin=(100,'K'), Tmax=(848.899,'K')), NASAPolynomial(coeffs=[1.93205,0.0558024,-2.44408e-05,4.49156e-09,-3.03791e-13,68999.2,33.6189], Tmin=(848.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C([CH2])C(9050)',
    structure = SMILES('[CH2]C[CH][CH]C([CH2])C'),
    E0 = (577.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,361.243,3063.57,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0257145,'amu*angstrom^2'), symmetry=1, barrier=(3.89672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257145,'amu*angstrom^2'), symmetry=1, barrier=(3.89672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257145,'amu*angstrom^2'), symmetry=1, barrier=(3.89672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257145,'amu*angstrom^2'), symmetry=1, barrier=(3.89672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257145,'amu*angstrom^2'), symmetry=1, barrier=(3.89672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257145,'amu*angstrom^2'), symmetry=1, barrier=(3.89672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25683,0.0491105,5.87504e-05,-2.44212e-07,2.25284e-10,69567.4,31.351], Tmin=(100,'K'), Tmax=(399.34,'K')), NASAPolynomial(coeffs=[3.64027,0.0530832,-2.31455e-05,4.31395e-09,-2.97176e-13,69314.7,24.1666], Tmin=(399.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C[CH]C[C]([CH2])C(5006)',
    structure = SMILES('[CH2]C[CH]C[C]([CH2])C'),
    E0 = (568.871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,485.572,2896.35,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0259788,'amu*angstrom^2'), symmetry=1, barrier=(4.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259788,'amu*angstrom^2'), symmetry=1, barrier=(4.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259788,'amu*angstrom^2'), symmetry=1, barrier=(4.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259788,'amu*angstrom^2'), symmetry=1, barrier=(4.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259788,'amu*angstrom^2'), symmetry=1, barrier=(4.03949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259788,'amu*angstrom^2'), symmetry=1, barrier=(4.03949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791942,0.0804866,-0.000105362,9.82149e-08,-3.6894e-11,68525.2,34.8391], Tmin=(100,'K'), Tmax=(851.22,'K')), NASAPolynomial(coeffs=[0.843925,0.0577459,-2.5646e-05,4.73465e-09,-3.20638e-13,69331.4,39.3841], Tmin=(851.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJ) + radical(Isobutyl)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3372.64,'J/mol'), sigma=(6.39812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.80 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.876276,0.0742873,-7.9907e-05,6.43255e-08,-2.21286e-11,69592,35.4413], Tmin=(100,'K'), Tmax=(870.302,'K')), NASAPolynomial(coeffs=[3.13779,0.05207,-2.12369e-05,3.77352e-09,-2.50661e-13,69646.1,27.4173], Tmin=(870.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (577.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1026.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (737.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (772.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (743.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1034.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1041.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1043.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1043.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1032.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1032.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (586.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (641.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (641.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (644.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (608.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (680.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (641.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (596.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (940.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (938.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (975.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (984.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (994.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (680.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (714.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (695.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (694.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (752.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (664.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (730.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (652.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (662.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (651.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (630.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['C3H6(27)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2][CH][CH]CC[CH2](655)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][CH][CH]CC[CH]C(3964)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][CH][CH]C[CH]CC(622)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]C([CH2])C([CH2])C(3893)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C[CH]C(556)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C([CH2])C(102)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][CH][C]CC([CH2])C(9042)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C][CH]CC([CH2])C(9043)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C(C)C[CH][CH][CH2](9044)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH][CH]CC([CH2])C(9045)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][CH]C1CC(C)C1(5009)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][C](C)CCC=C(5027)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2]C=C[CH]C(C)C(5017)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2]C=CC[C]([CH2])C(5011)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C=C[CH]C([CH2])C(5012)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C3H6(27)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00337229,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C3H6(T)(28)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C3H6(T)(28)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', '[CH2][CH][CH]C[CH][CH2](653)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH][CH]C[C]([CH2])C(9046)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C([CH2])C(9047)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.53274e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][CH][CH]C[C](C)C(9048)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C[CH]C([CH2])C(4096)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][CH][CH][CH]C(C)C(9049)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][CH]CC[C]([CH2])C(5005)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C[CH][CH]C([CH2])C(9050)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]CCC([CH2])[CH2](4139)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1855.84,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2]C[CH]C[C]([CH2])C(5006)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_1;Y_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2]C(C)[CH][CH][CH]C(3938)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(754000,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C[CH]CC([CH2])[CH2](608)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.69459e+07,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    products = ['[CH2][C](C)C[CH][CH]C(3937)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2][CH][CH]CC([CH2])C(3939)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(10873.3,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2013',
    isomers = [
        '[CH2][CH][CH]CC([CH2])C(3939)',
    ],
    reactants = [
        ('C3H6(27)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2013',
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


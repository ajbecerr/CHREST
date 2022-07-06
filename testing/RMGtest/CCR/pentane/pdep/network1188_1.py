species(
    label = '[CH2][CH]C(C)C[CH][O](1997)',
    structure = SMILES('[CH2][CH]C(C)C[CH][O]'),
    E0 = (451.752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,212.13,942.086,1324.36,1976.39],'cm^-1')),
        HinderedRotor(inertia=(0.112352,'amu*angstrom^2'), symmetry=1, barrier=(3.18431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112352,'amu*angstrom^2'), symmetry=1, barrier=(3.18431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112352,'amu*angstrom^2'), symmetry=1, barrier=(3.18431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112352,'amu*angstrom^2'), symmetry=1, barrier=(3.18431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112352,'amu*angstrom^2'), symmetry=1, barrier=(3.18431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751731,0.0793201,-0.000104001,9.31394e-08,-3.48471e-11,54442.6,32.9273], Tmin=(100,'K'), Tmax=(796.296,'K')), NASAPolynomial(coeffs=[3.83541,0.0502667,-2.37232e-05,4.53966e-09,-3.15353e-13,54381.5,21.4529], Tmin=(796.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = 'vinoxy(1351)',
    structure = SMILES('[CH2]C=O'),
    E0 = (6.53237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1127.76,1127.76,1127.77,1945],'cm^-1')),
        HinderedRotor(inertia=(0.00342653,'amu*angstrom^2'), symmetry=1, barrier=(3.09254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51396,0.0045106,2.69938e-05,-3.71579e-08,1.42465e-11,808.756,8.86322], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[6.49564,0.00770518,-2.52913e-06,4.69083e-10,-3.5128e-14,-479.661,-9.13857], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.53237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=CC(C)C[CH][O](1994)',
    structure = SMILES('C=CC(C)C[CH][O]'),
    E0 = (174.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,267.033,267.135,267.259,1282.08],'cm^-1')),
        HinderedRotor(inertia=(0.00236165,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00236229,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196822,'amu*angstrom^2'), symmetry=1, barrier=(9.96279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196473,'amu*angstrom^2'), symmetry=1, barrier=(9.96201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14439,0.0677673,-5.56557e-05,2.82003e-08,-6.54396e-12,21040.6,27.7015], Tmin=(100,'K'), Tmax=(969.587,'K')), NASAPolynomial(coeffs=[6.68749,0.0448996,-2.02785e-05,3.87591e-09,-2.72178e-13,19965.7,1.13034], Tmin=(969.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = '[CH2][CH]CC[CH][O](1487)',
    structure = SMILES('[CH2][CH]CC[CH][O]'),
    E0 = (477.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,1700.57,1701.02,1701.76],'cm^-1')),
        HinderedRotor(inertia=(0.13324,'amu*angstrom^2'), symmetry=1, barrier=(3.06346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13226,'amu*angstrom^2'), symmetry=1, barrier=(3.04093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132646,'amu*angstrom^2'), symmetry=1, barrier=(3.04979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132376,'amu*angstrom^2'), symmetry=1, barrier=(3.0436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3635.57,'J/mol'), sigma=(6.45581,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.87 K, Pc=30.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37842,0.0686508,-0.000108717,1.07913e-07,-4.11036e-11,57547.9,28.8918], Tmin=(100,'K'), Tmax=(853.943,'K')), NASAPolynomial(coeffs=[0.965186,0.0446752,-2.10876e-05,3.96924e-09,-2.704e-13,58563.2,36.3519], Tmin=(853.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH][O](1700)',
    structure = SMILES('[CH2]C([CH]C)C[CH][O]'),
    E0 = (451.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,248.404,890.758,1269.01,1987.23],'cm^-1')),
        HinderedRotor(inertia=(0.0940507,'amu*angstrom^2'), symmetry=1, barrier=(3.33577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0940507,'amu*angstrom^2'), symmetry=1, barrier=(3.33577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0940507,'amu*angstrom^2'), symmetry=1, barrier=(3.33577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0940507,'amu*angstrom^2'), symmetry=1, barrier=(3.33577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0940507,'amu*angstrom^2'), symmetry=1, barrier=(3.33577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.685187,0.0807927,-0.000107586,9.46849e-08,-3.41163e-11,54425.3,33.1627], Tmin=(100,'K'), Tmax=(830.633,'K')), NASAPolynomial(coeffs=[4.28037,0.0484164,-2.19173e-05,4.09419e-09,-2.79718e-13,54347.7,19.6132], Tmin=(830.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(C)[CH]C[CH][O](4364)',
    structure = SMILES('[CH2]C(C)[CH]C[CH][O]'),
    E0 = (448.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.685187,0.0807927,-0.000107586,9.46849e-08,-3.41163e-11,54022.7,33.1627], Tmin=(100,'K'), Tmax=(830.633,'K')), NASAPolynomial(coeffs=[4.28037,0.0484164,-2.19173e-05,4.09419e-09,-2.79718e-13,53945.1,19.6132], Tmin=(830.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])[O](2032)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])[O]'),
    E0 = (475.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,204.596,833.425,1690.42],'cm^-1')),
        HinderedRotor(inertia=(0.131802,'amu*angstrom^2'), symmetry=1, barrier=(3.28891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131802,'amu*angstrom^2'), symmetry=1, barrier=(3.28891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131802,'amu*angstrom^2'), symmetry=1, barrier=(3.28891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131802,'amu*angstrom^2'), symmetry=1, barrier=(3.28891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131802,'amu*angstrom^2'), symmetry=1, barrier=(3.28891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.824856,0.0704451,-5.62754e-05,2.453e-08,-4.52532e-12,57249.7,32.3517], Tmin=(100,'K'), Tmax=(1246.76,'K')), NASAPolynomial(coeffs=[11.2555,0.0369797,-1.6012e-05,3.00011e-09,-2.08094e-13,54648.8,-20.2706], Tmin=(1246.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(CJCO) + radical(RCCJ)"""),
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
    label = 'C[CH]C[CH][O](1560)',
    structure = SMILES('C[CH]C[CH][O]'),
    E0 = (296.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1857.88,1857.94],'cm^-1')),
        HinderedRotor(inertia=(0.102893,'amu*angstrom^2'), symmetry=1, barrier=(2.36571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103017,'amu*angstrom^2'), symmetry=1, barrier=(2.36856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102914,'amu*angstrom^2'), symmetry=1, barrier=(2.36619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23223,0.0293516,-1.06816e-05,9.14899e-10,6.61011e-14,35593.9,14.3951], Tmin=(100,'K'), Tmax=(2829.07,'K')), NASAPolynomial(coeffs=[40.9762,-0.00885552,2.28855e-06,-4.24163e-10,3.26746e-14,9303.22,-210.806], Tmin=(2829.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH)"""),
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
    label = '[CH2][C]C(C)C[CH][O](4365)',
    structure = SMILES('[CH2][C]C(C)C[CH][O]'),
    E0 = (705.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.500946,0.0836374,-0.000111987,9.38596e-08,-3.26413e-11,84962.7,31.1791], Tmin=(100,'K'), Tmax=(795.978,'K')), NASAPolynomial(coeffs=[6.92797,0.043085,-2.00111e-05,3.79697e-09,-2.62447e-13,84201,3.2816], Tmin=(795.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(705.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C(C)C[C][O](4366)',
    structure = SMILES('[CH2][CH]C(C)C[C][O]'),
    E0 = (732.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.832446,0.0762091,-9.67094e-05,8.23979e-08,-3.01877e-11,88200.3,31.8228], Tmin=(100,'K'), Tmax=(752.063,'K')), NASAPolynomial(coeffs=[5.39849,0.0452804,-2.17716e-05,4.22384e-09,-2.96687e-13,87701.4,12.3443], Tmin=(752.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(RCCJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH][CH]C(C)C[CH][O](4367)',
    structure = SMILES('[CH][CH]C(C)C[CH][O]'),
    E0 = (694.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645018,0.0818388,-0.000113705,1.01255e-07,-3.69802e-11,83668.7,32.5279], Tmin=(100,'K'), Tmax=(805.17,'K')), NASAPolynomial(coeffs=[5.04366,0.0462986,-2.19947e-05,4.20594e-09,-2.91507e-13,83404,15.0155], Tmin=(805.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1C(C)CC1[O](4368)',
    structure = SMILES('[CH2]C1C(C)CC1[O]'),
    E0 = (211.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36453,0.038569,4.99084e-05,-9.29641e-08,3.87475e-11,25514.3,25.1217], Tmin=(100,'K'), Tmax=(952.025,'K')), NASAPolynomial(coeffs=[15.5792,0.0273073,-8.70455e-06,1.55006e-09,-1.12859e-13,20611.5,-54.2915], Tmin=(952.025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C(C)CC[O](4369)',
    structure = SMILES('C=C[C](C)CC[O]'),
    E0 = (125.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48148,0.0563928,-2.94853e-05,7.23555e-09,-7.09696e-13,15219.3,25.169], Tmin=(100,'K'), Tmax=(2208.41,'K')), NASAPolynomial(coeffs=[14.9699,0.0319617,-1.28911e-05,2.22617e-09,-1.42616e-13,9261.75,-50.5915], Tmin=(2208.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Allyl_T)"""),
)

species(
    label = '[CH2]CC(C)C=C[O](4370)',
    structure = SMILES('[CH2]CC(C)C=C[O]'),
    E0 = (65.8622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655286,0.0580907,-2.33956e-06,-4.24513e-08,2.15634e-11,8055.74,29.3823], Tmin=(100,'K'), Tmax=(965.682,'K')), NASAPolynomial(coeffs=[17.5948,0.024984,-8.47872e-06,1.52679e-09,-1.09798e-13,3056.15,-60.6962], Tmin=(965.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.8622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C](C)CC=O(4183)',
    structure = SMILES('[CH2][CH][C](C)CC=O'),
    E0 = (304.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1487.97,2388.8],'cm^-1')),
        HinderedRotor(inertia=(0.11648,'amu*angstrom^2'), symmetry=1, barrier=(4.71688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00300222,'amu*angstrom^2'), symmetry=1, barrier=(4.71692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116483,'amu*angstrom^2'), symmetry=1, barrier=(4.71688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116477,'amu*angstrom^2'), symmetry=1, barrier=(4.71689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116483,'amu*angstrom^2'), symmetry=1, barrier=(4.71686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85395,0.0446575,-1.93424e-05,2.92602e-09,-1.11049e-13,36626.3,21.5782], Tmin=(100,'K'), Tmax=(2897.77,'K')), NASAPolynomial(coeffs=[55.9437,-0.0128829,3.0081e-06,-5.05532e-10,3.74369e-14,407.27,-295.55], Tmin=(2897.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C=O(4181)',
    structure = SMILES('[CH2][CH]C(C)C=C[O]'),
    E0 = (260.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,461.019,461.028,461.031],'cm^-1')),
        HinderedRotor(inertia=(0.000793129,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000793147,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110754,'amu*angstrom^2'), symmetry=1, barrier=(16.7042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110748,'amu*angstrom^2'), symmetry=1, barrier=(16.7042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835648,0.0569854,-1.40792e-05,-2.38866e-08,1.35378e-11,31444.5,31.2917], Tmin=(100,'K'), Tmax=(994.216,'K')), NASAPolynomial(coeffs=[15.6861,0.025602,-9.52352e-06,1.75347e-09,-1.24994e-13,27089.8,-47.3167], Tmin=(994.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][O](1556)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (367.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1667.93],'cm^-1')),
        HinderedRotor(inertia=(0.00517725,'amu*angstrom^2'), symmetry=1, barrier=(10.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22654,0.0212781,-3.59354e-05,3.91027e-08,-1.59281e-11,44278.5,12.2199], Tmin=(100,'K'), Tmax=(830.699,'K')), NASAPolynomial(coeffs=[2.17156,0.016018,-7.76586e-06,1.51127e-09,-1.05387e-13,44810.5,19.2612], Tmin=(830.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
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
    label = '[CH2]C=CC[CH][O](1848)',
    structure = SMILES('[CH2]C=CC[CH][O]'),
    E0 = (345.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,368.754,368.859,1940.62],'cm^-1')),
        HinderedRotor(inertia=(0.073507,'amu*angstrom^2'), symmetry=1, barrier=(7.09342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734922,'amu*angstrom^2'), symmetry=1, barrier=(7.09448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524237,'amu*angstrom^2'), symmetry=1, barrier=(50.5887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68117,0.0552208,-5.54951e-05,3.68237e-08,-1.10436e-11,41588.9,23.6157], Tmin=(100,'K'), Tmax=(778.818,'K')), NASAPolynomial(coeffs=[5.66422,0.0347641,-1.60961e-05,3.09862e-09,-2.18027e-13,40968.5,5.39543], Tmin=(778.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C[CH][O](1850)',
    structure = SMILES('[CH2][CH][CH]C[CH][O]'),
    E0 = (672.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180,1701.38,1701.66,1701.95],'cm^-1')),
        HinderedRotor(inertia=(0.195647,'amu*angstrom^2'), symmetry=1, barrier=(4.49832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218892,'amu*angstrom^2'), symmetry=1, barrier=(4.49774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218902,'amu*angstrom^2'), symmetry=1, barrier=(4.49819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195437,'amu*angstrom^2'), symmetry=1, barrier=(4.49349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63772,0.0690567,-0.000132374,1.43569e-07,-5.61224e-11,80920.8,30.5248], Tmin=(100,'K'), Tmax=(875.114,'K')), NASAPolynomial(coeffs=[-3.33206,0.0486816,-2.35884e-05,4.42807e-09,-2.98564e-13,83440.6,63.2656], Tmin=(875.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C](C)C[CH][O](4371)',
    structure = SMILES('[CH2][CH][C](C)C[CH][O]'),
    E0 = (637.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,453.99,1117.15,3442.72,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0477884,'amu*angstrom^2'), symmetry=1, barrier=(3.51222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477884,'amu*angstrom^2'), symmetry=1, barrier=(3.51222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477884,'amu*angstrom^2'), symmetry=1, barrier=(3.51222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477884,'amu*angstrom^2'), symmetry=1, barrier=(3.51222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477884,'amu*angstrom^2'), symmetry=1, barrier=(3.51222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7831,0.0874224,-0.000153399,1.58496e-07,-6.12707e-11,76734.2,33.9468], Tmin=(100,'K'), Tmax=(855.163,'K')), NASAPolynomial(coeffs=[-1.10166,0.0559408,-2.74946e-05,5.24033e-09,-3.58422e-13,78530.1,51.36], Tmin=(855.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(Cs_S) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][CH][O](4372)',
    structure = SMILES('[CH2][CH]C(C)[CH][CH][O]'),
    E0 = (651.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,551.497,1911.99,2185.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0654758,'amu*angstrom^2'), symmetry=1, barrier=(6.65298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0654758,'amu*angstrom^2'), symmetry=1, barrier=(6.65298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0654758,'amu*angstrom^2'), symmetry=1, barrier=(6.65298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0654758,'amu*angstrom^2'), symmetry=1, barrier=(6.65298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0654758,'amu*angstrom^2'), symmetry=1, barrier=(6.65298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939785,0.0741708,-9.55323e-05,8.37264e-08,-3.09645e-11,78479.5,34.7495], Tmin=(100,'K'), Tmax=(787.051,'K')), NASAPolynomial(coeffs=[4.4256,0.0459282,-2.16436e-05,4.14562e-09,-2.8843e-13,78256.8,20.8385], Tmin=(787.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCJCO) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH][O](4373)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH][O]'),
    E0 = (656.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,236.693,456.946,2042.85,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702035,0.0817735,-0.00011892,1.08066e-07,-3.90055e-11,79108.8,34.8881], Tmin=(100,'K'), Tmax=(846.966,'K')), NASAPolynomial(coeffs=[4.17012,0.0459529,-2.10488e-05,3.92751e-09,-2.67084e-13,79218.7,22.8492], Tmin=(846.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C](C)C[CH][O](4374)',
    structure = SMILES('[CH2]C[C](C)C[CH][O]'),
    E0 = (442.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59704,0.088565,-0.000141625,1.39598e-07,-5.29504e-11,53345.7,32.0595], Tmin=(100,'K'), Tmax=(849.195,'K')), NASAPolynomial(coeffs=[0.957344,0.0550748,-2.63101e-05,4.98118e-09,-3.40567e-13,54430.9,37.1297], Tmin=(849.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C[O](4375)',
    structure = SMILES('[CH2][CH]C(C)[CH]C[O]'),
    E0 = (471.357,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,189.209,472.954,2188.26,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0225382,'amu*angstrom^2'), symmetry=1, barrier=(2.90494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0225382,'amu*angstrom^2'), symmetry=1, barrier=(2.90494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0225382,'amu*angstrom^2'), symmetry=1, barrier=(2.90494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0225382,'amu*angstrom^2'), symmetry=1, barrier=(2.90494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0225382,'amu*angstrom^2'), symmetry=1, barrier=(2.90494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45262,0.0442551,5.82042e-05,-2.40434e-07,2.26823e-10,56735.9,30.1077], Tmin=(100,'K'), Tmax=(382.988,'K')), NASAPolynomial(coeffs=[3.35292,0.0499631,-2.3334e-05,4.51796e-09,-3.19095e-13,56556.2,25.1815], Tmin=(382.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.357,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(C)[CH][CH][O](4376)',
    structure = SMILES('[CH2]CC(C)[CH][CH][O]'),
    E0 = (457.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,221.923,833.05,1214.3,1775.38],'cm^-1')),
        HinderedRotor(inertia=(0.123273,'amu*angstrom^2'), symmetry=1, barrier=(3.48968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123273,'amu*angstrom^2'), symmetry=1, barrier=(3.48968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123273,'amu*angstrom^2'), symmetry=1, barrier=(3.48968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123273,'amu*angstrom^2'), symmetry=1, barrier=(3.48968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123273,'amu*angstrom^2'), symmetry=1, barrier=(3.48968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910219,0.0732342,-7.52654e-05,5.1788e-08,-1.60281e-11,55084.4,32.3152], Tmin=(100,'K'), Tmax=(758.877,'K')), NASAPolynomial(coeffs=[6.13299,0.0457045,-2.08485e-05,3.98189e-09,-2.78708e-13,54291.7,8.55954], Tmin=(758.877,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]CC([CH2])C[CH][O](1702)',
    structure = SMILES('[CH2]CC([CH2])C[CH][O]'),
    E0 = (462.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200.006,800.066,1199.91,1600.01],'cm^-1')),
        HinderedRotor(inertia=(0.156048,'amu*angstrom^2'), symmetry=1, barrier=(3.58842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156048,'amu*angstrom^2'), symmetry=1, barrier=(3.58842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156048,'amu*angstrom^2'), symmetry=1, barrier=(3.58842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156048,'amu*angstrom^2'), symmetry=1, barrier=(3.58842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156048,'amu*angstrom^2'), symmetry=1, barrier=(3.58842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.53007,0.0827358,-0.000106446,8.81565e-08,-3.02056e-11,55719.7,32.9512], Tmin=(100,'K'), Tmax=(829.817,'K')), NASAPolynomial(coeffs=[6.19,0.0451569,-1.9906e-05,3.67847e-09,-2.50084e-13,55134.9,8.83709], Tmin=(829.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][C](C)CC[O](4377)',
    structure = SMILES('[CH2][CH][C](C)CC[O]'),
    E0 = (456.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,187.635,1377.98,3174.23,3266.89],'cm^-1')),
        HinderedRotor(inertia=(0.0506916,'amu*angstrom^2'), symmetry=1, barrier=(1.77024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506916,'amu*angstrom^2'), symmetry=1, barrier=(1.77024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506916,'amu*angstrom^2'), symmetry=1, barrier=(1.77024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506916,'amu*angstrom^2'), symmetry=1, barrier=(1.77024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0506916,'amu*angstrom^2'), symmetry=1, barrier=(1.77024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19132,0.0451506,-1.71561e-05,1.49702e-09,1.17767e-13,54891.5,21.8783], Tmin=(100,'K'), Tmax=(2727.29,'K')), NASAPolynomial(coeffs=[53.98,-0.00754254,6.43645e-07,-1.20802e-10,1.55256e-14,20173.2,-282.138], Tmin=(2727.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][C](C)C[CH][O](1995)',
    structure = SMILES('C[CH][C](C)C[CH][O]'),
    E0 = (431.929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,180,1333.48,1988.77,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0923622,'amu*angstrom^2'), symmetry=1, barrier=(2.12359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0923622,'amu*angstrom^2'), symmetry=1, barrier=(2.12359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0923622,'amu*angstrom^2'), symmetry=1, barrier=(2.12359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0923622,'amu*angstrom^2'), symmetry=1, barrier=(2.12359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0923622,'amu*angstrom^2'), symmetry=1, barrier=(2.12359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755301,0.0865805,-0.000142598,1.4587e-07,-5.67311e-11,52051.1,32.26], Tmin=(100,'K'), Tmax=(847.431,'K')), NASAPolynomial(coeffs=[-0.957123,0.0583431,-2.83267e-05,5.39821e-09,-3.70312e-13,53645.5,47.9327], Tmin=(847.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(Cs_S) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][CH]O(4378)',
    structure = SMILES('[CH2][CH]C(C)[CH][CH]O'),
    E0 = (425.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.826139,0.0737683,-7.10424e-05,4.07207e-08,-1.00515e-11,51340.8,34.3106], Tmin=(100,'K'), Tmax=(951.557,'K')), NASAPolynomial(coeffs=[8.8755,0.0399323,-1.77055e-05,3.35314e-09,-2.34216e-13,49808.9,-4.1233], Tmin=(951.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(Cs_S) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])CC[O](4379)',
    structure = SMILES('[CH2][CH]C([CH2])CC[O]'),
    E0 = (476.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,417.822,443.829,2101.1,3694.3],'cm^-1')),
        HinderedRotor(inertia=(0.0242341,'amu*angstrom^2'), symmetry=1, barrier=(3.13072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242341,'amu*angstrom^2'), symmetry=1, barrier=(3.13072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242341,'amu*angstrom^2'), symmetry=1, barrier=(3.13072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242341,'amu*angstrom^2'), symmetry=1, barrier=(3.13072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242341,'amu*angstrom^2'), symmetry=1, barrier=(3.13072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985766,0.0728576,-8.71167e-05,7.51529e-08,-2.76886e-11,57416.3,34.3215], Tmin=(100,'K'), Tmax=(805.601,'K')), NASAPolynomial(coeffs=[3.46185,0.0493115,-2.23243e-05,4.19738e-09,-2.88961e-13,57382.5,25.1772], Tmin=(805.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C(C)[CH][CH][O](1996)',
    structure = SMILES('C[CH]C(C)[CH][CH][O]'),
    E0 = (446.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,203.661,793.213,913.899,3340.69],'cm^-1')),
        HinderedRotor(inertia=(0.0697783,'amu*angstrom^2'), symmetry=1, barrier=(1.91791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0697783,'amu*angstrom^2'), symmetry=1, barrier=(1.91791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0697783,'amu*angstrom^2'), symmetry=1, barrier=(1.91791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0697783,'amu*angstrom^2'), symmetry=1, barrier=(1.91791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0697783,'amu*angstrom^2'), symmetry=1, barrier=(1.91791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.962318,0.0726833,-8.2234e-05,6.7553e-08,-2.48024e-11,53794.3,32.8856], Tmin=(100,'K'), Tmax=(748.601,'K')), NASAPolynomial(coeffs=[4.3963,0.0486443,-2.26644e-05,4.34948e-09,-3.04224e-13,53439.6,18.3779], Tmin=(748.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][C](C)C[CH]O(4380)',
    structure = SMILES('[CH2][CH][C](C)C[CH]O'),
    E0 = (411.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44704,0.0898379,-0.00013983,1.31563e-07,-4.8277e-11,49605.1,34.2939], Tmin=(100,'K'), Tmax=(848.257,'K')), NASAPolynomial(coeffs=[3.54291,0.049606,-2.33584e-05,4.40062e-09,-3.00266e-13,50002,25.3035], Tmin=(848.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH]O(4381)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH]O'),
    E0 = (431.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.842673,0.077222,-7.27302e-05,2.1549e-08,1.10026e-11,51959.3,33.5999], Tmin=(100,'K'), Tmax=(586.415,'K')), NASAPolynomial(coeffs=[8.52837,0.0401657,-1.7255e-05,3.17349e-09,-2.16345e-13,50793.6,-1.63022], Tmin=(586.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(CCsJOH) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (451.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (451.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (896.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (594.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (646.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (632.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (909.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (917.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (917.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (944.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (906.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (460.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (515.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (515.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (524.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (479.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (526.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (515.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (493.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (818.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (808.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (848.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (863.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (868.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (617.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (608.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (573.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (580.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (573.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (600.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (613.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (552.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (582.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (575.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (511.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['vinoxy(1351)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['C=CC(C)C[CH][O](1994)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]CC[CH][O](1487)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2]C([CH]C)C[CH][O](1700)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2]C(C)[CH]C[CH][O](4364)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]C(C)C([CH2])[O](2032)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH2](502)', 'C[CH]C[CH][O](1560)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][O](1548)', '[CH2][CH]C([CH2])C(111)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]C(C)C[CH][O](4365)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH]C(C)C[C][O](4366)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH]C(C)C[CH][O](4367)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2]C1C(C)CC1[O](4368)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2]C=C(C)CC[O](4369)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2]CC(C)C=C[O](4370)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][CH][C](C)CC=O(4183)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2][CH]C(C)[CH]C=O(4181)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH][O](1556)', 'm1_allyl(186)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', '[CH2]C=CC[CH][O](1848)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.139979,'m^3/(mol*s)'), n=2.09962, Ea=(33.817,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['vinoxy(1351)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', '[CH2][CH][CH]C[CH][O](1850)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH][C](C)C[CH][O](4371)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][CH]C(C)[CH][CH][O](4372)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2]C[C](C)C[CH][O](4374)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C(C)[CH]C[O](4375)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CC(C)[CH][CH][O](4376)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CC([CH2])C[CH][O](1702)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][C](C)CC[O](4377)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['C[CH][C](C)C[CH][O](1995)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.081e+09,'s^-1'), n=0.812344, Ea=(148.932,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_Cs2] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2][CH]C(C)[CH][CH]O(4378)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.93361e+09,'s^-1'), n=1.0905, Ea=(161.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]C([CH2])CC[O](4379)'],
    products = ['[CH2][CH]C(C)C[CH][O](1997)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['C[CH]C(C)[CH][CH][O](1996)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2][CH][C](C)C[CH]O(4380)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(520772,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]C(C)C[CH][O](1997)'],
    products = ['[CH2][CH]C([CH2])C[CH]O(4381)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1188',
    isomers = [
        '[CH2][CH]C(C)C[CH][O](1997)',
    ],
    reactants = [
        ('vinoxy(1351)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1188',
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


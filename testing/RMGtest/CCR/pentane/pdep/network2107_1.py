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
    label = '[CH2][CH]C([CH2])CC=O(4184)',
    structure = SMILES('[CH2][CH]C([CH2])CC=O'),
    E0 = (324.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,180,1199.2],'cm^-1')),
        HinderedRotor(inertia=(0.176573,'amu*angstrom^2'), symmetry=1, barrier=(4.05977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176564,'amu*angstrom^2'), symmetry=1, barrier=(4.05955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00397902,'amu*angstrom^2'), symmetry=1, barrier=(4.06109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176523,'amu*angstrom^2'), symmetry=1, barrier=(4.0586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00397736,'amu*angstrom^2'), symmetry=1, barrier=(4.06151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7718,0.0570456,-2.14287e-05,-4.56858e-08,5.04468e-11,39102.1,30.0967], Tmin=(100,'K'), Tmax=(500.133,'K')), NASAPolynomial(coeffs=[5.11473,0.0433361,-1.93814e-05,3.66502e-09,-2.55097e-13,38604.8,14.6564], Tmin=(500.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH][O](4354)',
    structure = SMILES('[CH2][CH][CH]CC[CH][O]'),
    E0 = (648.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,476.327,863.886,2166.31,2892.04,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3778.62,'J/mol'), sigma=(6.76913,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=590.21 K, Pc=27.64 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944353,0.0844202,-0.000148429,1.54869e-07,-5.98592e-11,78085.5,35.2576], Tmin=(100,'K'), Tmax=(866.461,'K')), NASAPolynomial(coeffs=[-2.27782,0.0568642,-2.72688e-05,5.12854e-09,-3.47543e-13,80236.6,59.5321], Tmin=(866.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])[O](4207)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[O]'),
    E0 = (680.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,248.528,985.032,2220.78],'cm^-1')),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752264,0.073098,-7.18049e-05,4.06589e-08,-9.63206e-12,81917.4,34.4019], Tmin=(100,'K'), Tmax=(1005.11,'K')), NASAPolynomial(coeffs=[10.3935,0.0347274,-1.4539e-05,2.67408e-09,-1.83712e-13,79979.3,-12.1603], Tmin=(1005.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH][O](1558)',
    structure = SMILES('[CH2][CH]C[CH][O]'),
    E0 = (501.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,1713.42,1713.63],'cm^-1')),
        HinderedRotor(inertia=(0.153351,'amu*angstrom^2'), symmetry=1, barrier=(3.52584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00169109,'amu*angstrom^2'), symmetry=1, barrier=(3.52752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153281,'amu*angstrom^2'), symmetry=1, barrier=(3.52424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06734,0.053341,-9.28519e-05,9.68559e-08,-3.74628e-11,60383.4,24.1749], Tmin=(100,'K'), Tmax=(867.011,'K')), NASAPolynomial(coeffs=[-0.0668079,0.0364536,-1.73842e-05,3.26326e-09,-2.20958e-13,61758.3,39.9609], Tmin=(867.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
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
    label = '[CH2][C]C([CH2])C[CH][O](9297)',
    structure = SMILES('[CH2][C]C([CH2])C[CH][O]'),
    E0 = (910.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438455,0.0862649,-0.000127644,1.09984e-07,-3.74477e-11,109629,33.1844], Tmin=(100,'K'), Tmax=(860.077,'K')), NASAPolynomial(coeffs=[7.26422,0.0387672,-1.73338e-05,3.18401e-09,-2.14104e-13,109038,4.67005], Tmin=(860.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(910.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C[CH][O](9298)',
    structure = SMILES('[CH]C([CH][CH2])C[CH][O]'),
    E0 = (899.968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,205.712,920.246,1013.82,1182.79,1351.76,1601.45,1710.4],'cm^-1')),
        HinderedRotor(inertia=(0.145165,'amu*angstrom^2'), symmetry=1, barrier=(3.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145165,'amu*angstrom^2'), symmetry=1, barrier=(3.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145165,'amu*angstrom^2'), symmetry=1, barrier=(3.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145165,'amu*angstrom^2'), symmetry=1, barrier=(3.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145165,'amu*angstrom^2'), symmetry=1, barrier=(3.51319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655807,0.0828926,-0.000125292,1.14938e-07,-4.19704e-11,108352,34.2749], Tmin=(100,'K'), Tmax=(825.387,'K')), NASAPolynomial(coeffs=[4.96892,0.043772,-2.10886e-05,4.03014e-09,-2.78104e-13,108261,18.0536], Tmin=(825.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(899.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCsJOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[C][O](9299)',
    structure = SMILES('[CH2][CH]C([CH2])C[C][O]'),
    E0 = (937.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749322,0.0790861,-0.000113228,9.9505e-08,-3.52766e-11,112868,33.9014], Tmin=(100,'K'), Tmax=(826.447,'K')), NASAPolynomial(coeffs=[5.87223,0.0407157,-1.89462e-05,3.57488e-09,-2.45293e-13,112485,12.9678], Tmin=(826.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(937.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C[CH][O](9300)',
    structure = SMILES('[CH][CH]C([CH2])C[CH][O]'),
    E0 = (899.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,221.72,821.514,985.817,1150.12,1291.16,1594.9,1715.17],'cm^-1')),
        HinderedRotor(inertia=(0.132052,'amu*angstrom^2'), symmetry=1, barrier=(3.48897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132052,'amu*angstrom^2'), symmetry=1, barrier=(3.48897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132052,'amu*angstrom^2'), symmetry=1, barrier=(3.48897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132052,'amu*angstrom^2'), symmetry=1, barrier=(3.48897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132052,'amu*angstrom^2'), symmetry=1, barrier=(3.48897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600009,0.0842373,-0.000128443,1.15992e-07,-4.10949e-11,108335,34.4719], Tmin=(100,'K'), Tmax=(855.157,'K')), NASAPolynomial(coeffs=[5.34504,0.0420441,-1.93556e-05,3.60231e-09,-2.43957e-13,108254,16.5978], Tmin=(855.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(899.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC1C[CH][O](6570)',
    structure = SMILES('[CH2]C1CC1C[CH][O]'),
    E0 = (403.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.832445,0.0616348,-4.22325e-05,1.54018e-08,-2.3164e-12,48699.4,28.6695], Tmin=(100,'K'), Tmax=(1541.38,'K')), NASAPolynomial(coeffs=[13.2179,0.0294941,-1.0955e-05,1.87403e-09,-1.22336e-13,44881.2,-36.4421], Tmin=(1541.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1CC([O])C1[CH2](5528)',
    structure = SMILES('[CH2]C1CC([O])C1[CH2]'),
    E0 = (416.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39379,0.040108,3.80023e-05,-8.14092e-08,3.56403e-11,50177.1,26.799], Tmin=(100,'K'), Tmax=(931.218,'K')), NASAPolynomial(coeffs=[15.3623,0.0239487,-6.58837e-06,1.07106e-09,-7.57166e-14,45674.7,-49.8021], Tmin=(931.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CC([O])C1(5529)',
    structure = SMILES('[CH2][CH]C1CC([O])C1'),
    E0 = (414.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51773,0.0394184,2.88431e-05,-6.063e-08,2.43389e-11,50009.6,28.1078], Tmin=(100,'K'), Tmax=(1007.43,'K')), NASAPolynomial(coeffs=[12.879,0.0306433,-1.21923e-05,2.32672e-09,-1.68609e-13,45876.7,-35.9393], Tmin=(1007.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])CC=O(5534)',
    structure = SMILES('[CH2]C[C]([CH2])CC=O'),
    E0 = (315.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961988,0.076383,-0.000109865,1.0233e-07,-3.7678e-11,38033.2,31.4519], Tmin=(100,'K'), Tmax=(853.102,'K')), NASAPolynomial(coeffs=[2.56946,0.0475055,-2.15676e-05,4.00602e-09,-2.71593e-13,38535.5,28.5036], Tmin=(853.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(RCCJ) + radical(Isobutyl)"""),
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
    label = '[CH2]C=C([CH2])CC[O](9301)',
    structure = SMILES('[CH2]C=C([CH2])CC[O]'),
    E0 = (277.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19302,0.0632259,-4.45344e-05,1.67028e-08,-2.68721e-12,33447.7,26.8357], Tmin=(100,'K'), Tmax=(1377.86,'K')), NASAPolynomial(coeffs=[10.1186,0.0373145,-1.63261e-05,3.05438e-09,-2.10835e-13,30988,-19.086], Tmin=(1377.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])[CH]C=O(5535)',
    structure = SMILES('[CH2]CC([CH2])C=C[O]'),
    E0 = (270.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690581,0.0595632,-1.40415e-05,-3.1106e-08,1.85148e-11,32718.4,31.0378], Tmin=(100,'K'), Tmax=(938.115,'K')), NASAPolynomial(coeffs=[17.3253,0.0217135,-6.41292e-06,1.05961e-09,-7.36308e-14,28141.7,-55.9105], Tmin=(938.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][C]([CH2])CC=O(5532)',
    structure = SMILES('[CH2][CH][C]([CH2])CC=O'),
    E0 = (509.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1897.63],'cm^-1')),
        HinderedRotor(inertia=(0.142238,'amu*angstrom^2'), symmetry=1, barrier=(3.27033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00127978,'amu*angstrom^2'), symmetry=1, barrier=(3.27197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142231,'amu*angstrom^2'), symmetry=1, barrier=(3.27017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142081,'amu*angstrom^2'), symmetry=1, barrier=(3.26673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00127868,'amu*angstrom^2'), symmetry=1, barrier=(3.26639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14813,0.0752401,-0.000121643,1.21243e-07,-4.6011e-11,61421.7,33.3389], Tmin=(100,'K'), Tmax=(860.165,'K')), NASAPolynomial(coeffs=[0.507294,0.048377,-2.27554e-05,4.26594e-09,-2.89512e-13,62636,42.7515], Tmin=(860.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH]C=O(5533)',
    structure = SMILES('[CH2][CH]C([CH2])C=C[O]'),
    E0 = (465.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,726.327,3368.12],'cm^-1')),
        HinderedRotor(inertia=(0.0081838,'amu*angstrom^2'), symmetry=1, barrier=(3.08188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8863,'amu*angstrom^2'), symmetry=1, barrier=(20.3778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0547218,'amu*angstrom^2'), symmetry=1, barrier=(20.377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219919,'amu*angstrom^2'), symmetry=1, barrier=(82.3948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.872014,0.0584622,-2.5888e-05,-1.2245e-08,1.02813e-11,56107,32.9424], Tmin=(100,'K'), Tmax=(960.31,'K')), NASAPolynomial(coeffs=[15.3348,0.0224657,-7.53297e-06,1.30371e-09,-9.02495e-14,52211.3,-42.0672], Tmin=(960.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = 'C=CC[CH][O](1525)',
    structure = SMILES('C=CC[CH][O]'),
    E0 = (229.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,245.685,245.711,1409.83],'cm^-1')),
        HinderedRotor(inertia=(0.128793,'amu*angstrom^2'), symmetry=1, barrier=(5.51798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128739,'amu*angstrom^2'), symmetry=1, barrier=(5.51754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16344,0.0442977,-5.21517e-05,4.48084e-08,-1.67955e-11,27683.3,19.7843], Tmin=(100,'K'), Tmax=(774.377,'K')), NASAPolynomial(coeffs=[3.80673,0.0302191,-1.40523e-05,2.68592e-09,-1.87066e-13,27596.4,13.3588], Tmin=(774.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = '[CH2][CH][C]([CH2])C[CH][O](9302)',
    structure = SMILES('[CH2][CH][C]([CH2])C[CH][O]'),
    E0 = (842.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,194.99,204.443,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0099366,'amu*angstrom^2'), symmetry=1, barrier=(7.01745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0099366,'amu*angstrom^2'), symmetry=1, barrier=(7.01745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0099366,'amu*angstrom^2'), symmetry=1, barrier=(7.01745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0099366,'amu*angstrom^2'), symmetry=1, barrier=(7.01745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0099366,'amu*angstrom^2'), symmetry=1, barrier=(7.01745,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773362,0.0893917,-0.000166605,1.71286e-07,-6.46231e-11,101399,35.7654], Tmin=(100,'K'), Tmax=(877.018,'K')), NASAPolynomial(coeffs=[-0.976732,0.0519978,-2.50401e-05,4.68117e-09,-3.14615e-13,103451,53.9278], Tmin=(877.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(842.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(Cs_S) + radical(CCsJOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][CH][O](9303)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH][O]'),
    E0 = (856.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,235.884,2407.44,2826.16,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0489722,'amu*angstrom^2'), symmetry=1, barrier=(5.76455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0489722,'amu*angstrom^2'), symmetry=1, barrier=(5.76455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0489722,'amu*angstrom^2'), symmetry=1, barrier=(5.76455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0489722,'amu*angstrom^2'), symmetry=1, barrier=(5.76455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0489722,'amu*angstrom^2'), symmetry=1, barrier=(5.76455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882772,0.0767179,-0.000110812,9.91639e-08,-3.53564e-11,103146,36.7361], Tmin=(100,'K'), Tmax=(846.803,'K')), NASAPolynomial(coeffs=[4.78482,0.0415703,-1.89427e-05,3.52701e-09,-2.39613e-13,103084,22.0985], Tmin=(846.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(856.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCJCO) + radical(Isobutyl) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C[C]([CH2])C[CH][O](9304)',
    structure = SMILES('[CH2]C[C]([CH2])C[CH][O]'),
    E0 = (647.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580101,0.0906237,-0.000155163,1.52841e-07,-5.6503e-11,78010.6,33.9036], Tmin=(100,'K'), Tmax=(875.557,'K')), NASAPolynomial(coeffs=[1.10941,0.051084,-2.38273e-05,4.41521e-09,-2.96188e-13,79340.7,39.5458], Tmin=(875.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(CCsJOH) + radical(RCCJ) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C([CH2])[CH]C[O](5522)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]C[O]'),
    E0 = (676.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1065.48,1965.12,3476.42,4000],'cm^-1')),
        HinderedRotor(inertia=(0.106922,'amu*angstrom^2'), symmetry=1, barrier=(5.1716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106922,'amu*angstrom^2'), symmetry=1, barrier=(5.1716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106922,'amu*angstrom^2'), symmetry=1, barrier=(5.1716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106922,'amu*angstrom^2'), symmetry=1, barrier=(5.1716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106922,'amu*angstrom^2'), symmetry=1, barrier=(5.1716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17399,0.0677053,-7.86295e-05,6.56988e-08,-2.37767e-11,81453.2,36.1432], Tmin=(100,'K'), Tmax=(796.624,'K')), NASAPolynomial(coeffs=[4.0551,0.0449675,-2.02414e-05,3.80253e-09,-2.61971e-13,81256.6,24.5458], Tmin=(796.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[CH][CH][O](9305)',
    structure = SMILES('[CH2]CC([CH2])[CH][CH][O]'),
    E0 = (662.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,232.1,1033.45,1506.25,2142.92],'cm^-1')),
        HinderedRotor(inertia=(0.0651954,'amu*angstrom^2'), symmetry=1, barrier=(2.50147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0651954,'amu*angstrom^2'), symmetry=1, barrier=(2.50147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0651954,'amu*angstrom^2'), symmetry=1, barrier=(2.50147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0651954,'amu*angstrom^2'), symmetry=1, barrier=(2.50147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0651954,'amu*angstrom^2'), symmetry=1, barrier=(2.50147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.714235,0.0776358,-9.81628e-05,7.89953e-08,-2.64302e-11,79756.8,34.7871], Tmin=(100,'K'), Tmax=(825.51,'K')), NASAPolynomial(coeffs=[6.79647,0.040789,-1.78087e-05,3.28012e-09,-2.22796e-13,79003.9,8.1322], Tmin=(825.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(Isobutyl) + radical(RCCJ) + radical(CCsJOH)"""),
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
    label = '[CH2][CH][C]([CH2])CC[O](5523)',
    structure = SMILES('[CH2][CH][C]([CH2])CC[O]'),
    E0 = (661.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2262.13,2374.56,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.025753,'amu*angstrom^2'), symmetry=1, barrier=(12.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025753,'amu*angstrom^2'), symmetry=1, barrier=(12.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025753,'amu*angstrom^2'), symmetry=1, barrier=(12.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025753,'amu*angstrom^2'), symmetry=1, barrier=(12.3182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025753,'amu*angstrom^2'), symmetry=1, barrier=(12.3182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02071,0.0809251,-0.000136447,1.40554e-07,-5.42216e-11,79707.8,35.3278], Tmin=(100,'K'), Tmax=(863.489,'K')), NASAPolynomial(coeffs=[-1.52418,0.0550718,-2.61468e-05,4.91032e-09,-3.33063e-13,81550.6,55.3578], Tmin=(863.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C[CH][O](4389)',
    structure = SMILES('[CH2][C]([CH]C)C[CH][O]'),
    E0 = (637.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,979.774,2304.5,2826.56,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0762111,'amu*angstrom^2'), symmetry=1, barrier=(7.09362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0762111,'amu*angstrom^2'), symmetry=1, barrier=(7.09362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0762111,'amu*angstrom^2'), symmetry=1, barrier=(7.09362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0762111,'amu*angstrom^2'), symmetry=1, barrier=(7.09362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0762111,'amu*angstrom^2'), symmetry=1, barrier=(7.09362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741465,0.0885968,-0.000155954,1.58818e-07,-6.01241e-11,76715.8,34.0933], Tmin=(100,'K'), Tmax=(871.375,'K')), NASAPolynomial(coeffs=[-0.802577,0.0543485,-2.58418e-05,4.83178e-09,-3.25897e-13,78554.2,50.3347], Tmin=(871.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Tertalkyl) + radical(Cs_S) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][CH]O(9306)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH]O'),
    E0 = (631.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598617,0.0784587,-9.45643e-05,6.82452e-08,-2.0405e-11,76014.7,36.901], Tmin=(100,'K'), Tmax=(810.3,'K')), NASAPolynomial(coeffs=[9.30871,0.035453,-1.4937e-05,2.71906e-09,-1.8415e-13,74603.4,-3.28624], Tmin=(810.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][CH][O](4390)',
    structure = SMILES('[CH2]C([CH]C)[CH][CH][O]'),
    E0 = (651.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,224.01,992.843,1615.35,3075.99],'cm^-1')),
        HinderedRotor(inertia=(0.0601551,'amu*angstrom^2'), symmetry=1, barrier=(2.36532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601551,'amu*angstrom^2'), symmetry=1, barrier=(2.36532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601551,'amu*angstrom^2'), symmetry=1, barrier=(2.36532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601551,'amu*angstrom^2'), symmetry=1, barrier=(2.36532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601551,'amu*angstrom^2'), symmetry=1, barrier=(2.36532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867928,0.0757113,-9.93775e-05,8.5636e-08,-3.03969e-11,78462.4,35.0036], Tmin=(100,'K'), Tmax=(827.647,'K')), NASAPolynomial(coeffs=[4.88973,0.0440433,-1.98169e-05,3.69508e-09,-2.52365e-13,78215.6,18.8922], Tmin=(827.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCJCO) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C[CH]O(5524)',
    structure = SMILES('[CH2][CH][C]([CH2])C[CH]O'),
    E0 = (616.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,236.038,2057.23,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0542525,'amu*angstrom^2'), symmetry=1, barrier=(2.25739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542525,'amu*angstrom^2'), symmetry=1, barrier=(2.25739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542525,'amu*angstrom^2'), symmetry=1, barrier=(2.25739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542525,'amu*angstrom^2'), symmetry=1, barrier=(2.25739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542525,'amu*angstrom^2'), symmetry=1, barrier=(2.25739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542525,'amu*angstrom^2'), symmetry=1, barrier=(2.25739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.426,0.0919506,-0.000153586,1.4514e-07,-5.19993e-11,74270.1,36.1524], Tmin=(100,'K'), Tmax=(879.279,'K')), NASAPolynomial(coeffs=[3.69933,0.0456072,-2.08706e-05,3.83343e-09,-2.55782e-13,74910.3,27.6955], Tmin=(879.279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(616.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(CCsJOH) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (656.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (656.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (814.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (837.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1114.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1110.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1122.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1122.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1111.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1149.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1111.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (662.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (665.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (665.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (679.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (679.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (720.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (720.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (720.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (729.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (685.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (766.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (787.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (679.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (698.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1023.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1054.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1068.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (815.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (798.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (813.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (778.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (774.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (778.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (814.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (812.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (787.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (780.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['vinoxy(1351)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH]C([CH2])CC=O(4184)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C[CH][O](1558)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C[CH][O](1850)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][O](1548)', '[CH2][CH]C([CH2])[CH2](499)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C]C([CH2])C[CH][O](9297)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([CH][CH2])C[CH][O](9298)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[C][O](9299)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH]C([CH2])C[CH][O](9300)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]C1CC1C[CH][O](6570)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]C1CC([O])C1[CH2](5528)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH]C1CC([O])C1(5529)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]C[C]([CH2])CC=O(5534)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH][C](C)CC=O(4183)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]C=C([CH2])CC[O](9301)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]CC([CH2])[CH]C=O(5535)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH]C(C)[CH]C=O(4181)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])CC=O(5532)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[CH]C=O(5533)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(T)(20)', '[CH2]C=CC[CH][O](1848)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH][CH2](502)', 'C=CC[CH][O](1525)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['vinoxy(1351)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C[CH][O](9302)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[CH][CH][O](9303)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]C[C]([CH2])C[CH][O](9304)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH][C](C)C[CH][O](4371)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]C([CH2])[CH]C[O](5522)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]CC([CH2])[CH][CH][O](9305)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH]C(C)[CH][CH][O](4372)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH][C]([CH2])CC[O](5523)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][C]([CH]C)C[CH][O](4389)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH]C([CH2])[CH][CH]O(9306)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.30814e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2]C([CH]C)[CH][CH][O](4390)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH][C]([CH2])C[CH]O(5524)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(520772,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2107',
    isomers = [
        '[CH2][CH]C([CH2])C[CH][O](4373)',
    ],
    reactants = [
        ('vinoxy(1351)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2107',
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


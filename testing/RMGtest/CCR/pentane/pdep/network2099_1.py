species(
    label = '[CH2][CH][CH]CC([CH2])[O](4458)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[O]'),
    E0 = (671.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,686.59,2061.2,2278.4,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852795,0.0775059,-0.000107974,9.70252e-08,-3.50878e-11,80900.1,35.2748], Tmin=(100,'K'), Tmax=(840.862,'K')), NASAPolynomial(coeffs=[4.01348,0.0455248,-2.06945e-05,3.86014e-09,-2.62948e-13,80967.6,24.1364], Tmin=(840.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
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
    label = '[CH2]C([CH2])[O](1366)',
    structure = SMILES('[CH2]C([CH2])[O]'),
    E0 = (362.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,450.564],'cm^-1')),
        HinderedRotor(inertia=(0.000831827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.050104,'amu*angstrom^2'), symmetry=1, barrier=(7.38455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17517,0.0359905,-3.15384e-05,1.43889e-08,-2.59785e-12,43692.5,16.9615], Tmin=(100,'K'), Tmax=(1341.41,'K')), NASAPolynomial(coeffs=[10.4304,0.0113739,-4.01123e-06,7.08062e-10,-4.81153e-14,41477.8,-25.2898], Tmin=(1341.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CC(C)OJ) + radical(CJCO) + radical(CJCO)"""),
)

species(
    label = 'O(4)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
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
    label = '[CH2][CH][C]CC([CH2])[O](9315)',
    structure = SMILES('[CH2][CH][C]CC([CH2])[O]'),
    E0 = (925.523,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81903,0.0762403,-8.64556e-05,5.14094e-08,-9.44319e-12,111423,32.7749], Tmin=(100,'K'), Tmax=(641.51,'K')), NASAPolynomial(coeffs=[8.89739,0.0360246,-1.61672e-05,3.04096e-09,-2.1026e-13,110178,-4.24095], Tmin=(641.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(925.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJC) + radical(CJCO) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC([CH2])[O](9316)',
    structure = SMILES('[CH2][C][CH]CC([CH2])[O]'),
    E0 = (925.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60653,0.0802528,-0.000107462,8.70094e-08,-2.90537e-11,111433,33.3314], Tmin=(100,'K'), Tmax=(800.107,'K')), NASAPolynomial(coeffs=[8.08223,0.0374541,-1.70542e-05,3.20494e-09,-2.20294e-13,110410,0.0178419], Tmin=(800.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(925.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(CJCO) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])C[CH][CH][CH2](9018)',
    structure = SMILES('[CH]C([O])C[CH][CH][CH2]'),
    E0 = (908.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,216.35,1045.48,1271.49,1332.58,1522.69,1813.71,2048.82],'cm^-1')),
        HinderedRotor(inertia=(0.0831084,'amu*angstrom^2'), symmetry=1, barrier=(2.92489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0831084,'amu*angstrom^2'), symmetry=1, barrier=(2.92489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0831084,'amu*angstrom^2'), symmetry=1, barrier=(2.92489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0831084,'amu*angstrom^2'), symmetry=1, barrier=(2.92489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0831084,'amu*angstrom^2'), symmetry=1, barrier=(2.92489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950666,0.0753968,-0.000107829,9.76543e-08,-3.53507e-11,109356,35.2937], Tmin=(100,'K'), Tmax=(842.231,'K')), NASAPolynomial(coeffs=[4.09497,0.043228,-1.98401e-05,3.70959e-09,-2.52793e-13,109438,24.2932], Tmin=(842.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(908.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC([CH2])[O](9317)',
    structure = SMILES('[CH][CH][CH]CC([CH2])[O]'),
    E0 = (914.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,211.481,1045.19,1054.62,1256.48,1402.42,1706.44,1866.89],'cm^-1')),
        HinderedRotor(inertia=(0.114846,'amu*angstrom^2'), symmetry=1, barrier=(3.22053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114846,'amu*angstrom^2'), symmetry=1, barrier=(3.22053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114846,'amu*angstrom^2'), symmetry=1, barrier=(3.22053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114846,'amu*angstrom^2'), symmetry=1, barrier=(3.22053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114846,'amu*angstrom^2'), symmetry=1, barrier=(3.22053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749783,0.0799816,-0.00011754,1.05005e-07,-3.71981e-11,110126,34.8622], Tmin=(100,'K'), Tmax=(849.9,'K')), NASAPolynomial(coeffs=[5.19354,0.0416069,-1.89959e-05,3.53364e-09,-2.39712e-13,110001,17.8562], Tmin=(849.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(914.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][CH]CC1CO1(9318)',
    structure = SMILES('[CH2][CH][CH]CC1CO1'),
    E0 = (419.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23906,0.0623757,-6.04531e-05,4.16397e-08,-1.21382e-11,50537.5,30.4858], Tmin=(100,'K'), Tmax=(981.403,'K')), NASAPolynomial(coeffs=[5.11123,0.0391955,-1.37166e-05,2.21059e-09,-1.37459e-13,50133.8,13.6926], Tmin=(981.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])O1(5565)',
    structure = SMILES('[CH2][CH]C1CC([CH2])O1'),
    E0 = (422.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05909,0.0536278,-2.12268e-06,-4.32976e-08,2.43253e-11,50957,27.4182], Tmin=(100,'K'), Tmax=(877.909,'K')), NASAPolynomial(coeffs=[14.7512,0.0233221,-5.15316e-06,6.25977e-10,-3.54616e-14,47316.7,-43.8959], Tmin=(877.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJCO) + radical(CJC(C)OC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]CC(=C)O(4326)',
    structure = SMILES('[CH2][CH][CH]CC(=C)O'),
    E0 = (315.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,203.845,522.076,2865.65],'cm^-1')),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513497,'amu*angstrom^2'), symmetry=1, barrier=(1.45703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03723,0.0676371,-6.46233e-05,3.7456e-08,-9.31459e-12,38051.7,32.3541], Tmin=(100,'K'), Tmax=(949.571,'K')), NASAPolynomial(coeffs=[8.4528,0.0363987,-1.52762e-05,2.81001e-09,-1.9288e-13,36643.4,-3.03783], Tmin=(949.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)=O(4331)',
    structure = SMILES('[CH2][CH][CH]CC(C)=O'),
    E0 = (293.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3139,0.0642699,-7.04167e-05,5.80813e-08,-2.14531e-11,35432.6,32.134], Tmin=(100,'K'), Tmax=(769.892,'K')), NASAPolynomial(coeffs=[3.70612,0.0452824,-2.06448e-05,3.91782e-09,-2.72106e-13,35258.6,22.4809], Tmin=(769.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])O(5620)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])O'),
    E0 = (231.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0720449,0.0786114,-6.94006e-05,3.11102e-08,-5.51262e-12,27963.5,29.2252], Tmin=(100,'K'), Tmax=(1366.72,'K')), NASAPolynomial(coeffs=[18.8038,0.0233672,-8.76926e-06,1.53519e-09,-1.02792e-13,22803.9,-67.7372], Tmin=(1366.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)[O](5621)',
    structure = SMILES('[CH2]C=C[CH]C(C)[O]'),
    E0 = (249.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,421.193,421.653,422.395],'cm^-1')),
        HinderedRotor(inertia=(0.203033,'amu*angstrom^2'), symmetry=1, barrier=(25.5492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704557,'amu*angstrom^2'), symmetry=1, barrier=(89.151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202815,'amu*angstrom^2'), symmetry=1, barrier=(25.56,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201783,'amu*angstrom^2'), symmetry=1, barrier=(25.5538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420245,0.0674528,-3.73873e-05,-5.16328e-10,5.06711e-12,30204.1,26.6255], Tmin=(100,'K'), Tmax=(1048.99,'K')), NASAPolynomial(coeffs=[16.3814,0.0278217,-1.10777e-05,2.05836e-09,-1.45074e-13,25687.3,-56.7091], Tmin=(1048.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CCC(=C)[O](4322)',
    structure = SMILES('[CH2][CH]CCC(=C)[O]'),
    E0 = (258.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,283.133,283.138,1104.16,2570.2],'cm^-1')),
        HinderedRotor(inertia=(0.17237,'amu*angstrom^2'), symmetry=1, barrier=(9.80502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0021029,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210292,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.075,0.0641934,-5.22537e-05,2.4342e-08,-4.83332e-12,31238.6,31.147], Tmin=(100,'K'), Tmax=(1169.1,'K')), NASAPolynomial(coeffs=[9.778,0.0344163,-1.40481e-05,2.55539e-09,-1.74434e-13,29203.7,-12.1996], Tmin=(1169.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])[O](5434)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])[O]'),
    E0 = (461.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,459.869,459.869,459.87],'cm^-1')),
        HinderedRotor(inertia=(0.598018,'amu*angstrom^2'), symmetry=1, barrier=(89.7449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179001,'amu*angstrom^2'), symmetry=1, barrier=(26.8625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178999,'amu*angstrom^2'), symmetry=1, barrier=(26.8625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179,'amu*angstrom^2'), symmetry=1, barrier=(26.8625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234313,0.0730297,-5.84454e-05,2.02624e-08,-1.74389e-12,55657.7,27.9272], Tmin=(100,'K'), Tmax=(1090.56,'K')), NASAPolynomial(coeffs=[17.3969,0.0237107,-9.35774e-06,1.71542e-09,-1.1942e-13,51103.8,-60.0762], Tmin=(1090.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]CC(=C)[O](4319)',
    structure = SMILES('[CH2][CH][CH]CC(=C)[O]'),
    E0 = (453.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1691.94,1692.09,1692.41,1692.43],'cm^-1')),
        HinderedRotor(inertia=(0.0775941,'amu*angstrom^2'), symmetry=1, barrier=(5.04212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0773371,'amu*angstrom^2'), symmetry=1, barrier=(5.03386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.078018,'amu*angstrom^2'), symmetry=1, barrier=(5.06108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0779554,'amu*angstrom^2'), symmetry=1, barrier=(5.04509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2439,0.0656708,-7.98273e-05,6.56523e-08,-2.26796e-11,54615.6,33.1053], Tmin=(100,'K'), Tmax=(830.908,'K')), NASAPolynomial(coeffs=[4.97977,0.0392941,-1.70605e-05,3.13676e-09,-2.12872e-13,54284.4,17.5173], Tmin=(830.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH][CH]C([CH2])[O](9139)',
    structure = SMILES('[CH2][CH][CH][CH]C([CH2])[O]'),
    E0 = (871.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1976.85,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.037357,'amu*angstrom^2'), symmetry=1, barrier=(15.9284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037357,'amu*angstrom^2'), symmetry=1, barrier=(15.9284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037357,'amu*angstrom^2'), symmetry=1, barrier=(15.9284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037357,'amu*angstrom^2'), symmetry=1, barrier=(15.9284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037357,'amu*angstrom^2'), symmetry=1, barrier=(15.9284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0341,0.072443,-9.98372e-05,8.80813e-08,-3.14185e-11,104937,37.1208], Tmin=(100,'K'), Tmax=(839.704,'K')), NASAPolynomial(coeffs=[4.62679,0.0411447,-1.85898e-05,3.46e-09,-2.35508e-13,104834,23.3934], Tmin=(839.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(871.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C]([CH2])[O](9319)',
    structure = SMILES('[CH2][CH][CH]C[C]([CH2])[O]'),
    E0 = (848.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,199.027,220.449,3403.69,4000],'cm^-1')),
        HinderedRotor(inertia=(0.012835,'amu*angstrom^2'), symmetry=1, barrier=(6.56263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012835,'amu*angstrom^2'), symmetry=1, barrier=(6.56263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012835,'amu*angstrom^2'), symmetry=1, barrier=(6.56263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012835,'amu*angstrom^2'), symmetry=1, barrier=(6.56263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012835,'amu*angstrom^2'), symmetry=1, barrier=(6.56263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.66998,0.0865883,-0.000146273,1.40866e-07,-5.13506e-11,102145,35.3545], Tmin=(100,'K'), Tmax=(872.384,'K')), NASAPolynomial(coeffs=[2.81653,0.0452888,-2.11736e-05,3.93168e-09,-2.64266e-13,102968,32.1524], Tmin=(872.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(848.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(C2CsJOH) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])[O](5614)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])[O]'),
    E0 = (677.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,228.938,1130.48,1456.53,2126.86],'cm^-1')),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0803259,'amu*angstrom^2'), symmetry=1, barrier=(2.94639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00051,0.0692042,-6.53765e-05,3.68934e-08,-8.9716e-12,81554.8,34.6887], Tmin=(100,'K'), Tmax=(965.242,'K')), NASAPolynomial(coeffs=[8.60266,0.0377007,-1.64197e-05,3.08033e-09,-2.13975e-13,80087.2,-1.71839], Tmin=(965.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C]([CH2])O(9320)',
    structure = SMILES('[CH2][CH][CH]C[C]([CH2])O'),
    E0 = (618.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509362,0.0905514,-0.000152043,1.45666e-07,-5.28223e-11,74444.7,36.1231], Tmin=(100,'K'), Tmax=(875.904,'K')), NASAPolynomial(coeffs=[2.89845,0.0470528,-2.17425e-05,4.01515e-09,-2.68915e-13,75276.3,32.0498], Tmin=(875.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)[O](9321)',
    structure = SMILES('[CH2][CH][CH]C[C](C)[O]'),
    E0 = (636.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.840219,0.0811923,-0.000125845,1.20945e-07,-4.49509e-11,76692.4,34.1096], Tmin=(100,'K'), Tmax=(858.481,'K')), NASAPolynomial(coeffs=[1.89118,0.0492622,-2.28204e-05,4.25847e-09,-2.88647e-13,77508.1,35.0015], Tmin=(858.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C]([CH2])[O](5615)',
    structure = SMILES('[CH2][CH]CC[C]([CH2])[O]'),
    E0 = (653.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,232.124,807.996,1204.69,1799.26],'cm^-1')),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123174,'amu*angstrom^2'), symmetry=1, barrier=(3.35367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430242,0.0859354,-0.000121676,1.03884e-07,-3.57202e-11,78771.7,33.6524], Tmin=(100,'K'), Tmax=(837.086,'K')), NASAPolynomial(coeffs=[7.05161,0.0413927,-1.87383e-05,3.48865e-09,-2.37434e-13,78115.2,5.58572], Tmin=(837.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C([CH2])O(9322)',
    structure = SMILES('[CH2][CH][CH][CH]C([CH2])O'),
    E0 = (641.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.870571,0.0764405,-0.000105724,9.30185e-08,-3.29362e-11,77236.8,37.8998], Tmin=(100,'K'), Tmax=(847.677,'K')), NASAPolynomial(coeffs=[4.72713,0.0428763,-1.91397e-05,3.53889e-09,-2.39773e-13,77135,23.1878], Tmin=(847.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)[O](9323)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)[O]'),
    E0 = (660.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22907,0.0667403,-7.82808e-05,6.66565e-08,-2.43837e-11,79483.3,35.7882], Tmin=(100,'K'), Tmax=(806.029,'K')), NASAPolynomial(coeffs=[3.593,0.0453103,-2.0351e-05,3.81436e-09,-2.62215e-13,79417.3,26.8478], Tmin=(806.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C([CH2])[O](9324)',
    structure = SMILES('[CH2]C[CH][CH]C([CH2])[O]'),
    E0 = (677.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,275.04,782.909,1731.27,3036.06],'cm^-1')),
        HinderedRotor(inertia=(0.0504303,'amu*angstrom^2'), symmetry=1, barrier=(2.33681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0504303,'amu*angstrom^2'), symmetry=1, barrier=(2.33681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0504303,'amu*angstrom^2'), symmetry=1, barrier=(2.33681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0504303,'amu*angstrom^2'), symmetry=1, barrier=(2.33681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0504303,'amu*angstrom^2'), symmetry=1, barrier=(2.33681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02933,0.0697143,-7.0389e-05,4.46532e-08,-1.23917e-11,81554.1,34.4176], Tmin=(100,'K'), Tmax=(848.939,'K')), NASAPolynomial(coeffs=[7.34342,0.0399636,-1.78217e-05,3.37217e-09,-2.35012e-13,80482,4.98982], Tmin=(848.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJCC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[C]([CH2])[O](5616)',
    structure = SMILES('[CH2]C[CH]C[C]([CH2])[O]'),
    E0 = (653.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,185.833,675.106,843.071,3300.92],'cm^-1')),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947527,'amu*angstrom^2'), symmetry=1, barrier=(2.25048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467153,0.0864175,-0.000126902,1.1233e-07,-3.95982e-11,78770.5,33.3478], Tmin=(100,'K'), Tmax=(841.968,'K')), NASAPolynomial(coeffs=[5.97125,0.0433223,-1.99352e-05,3.72974e-09,-2.54111e-13,78444.3,11.3076], Tmin=(841.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])[CH][CH][CH]C(4312)',
    structure = SMILES('[CH2]C([O])[CH][CH][CH]C'),
    E0 = (666.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,217.072,360.097,2368.32,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0165365,'amu*angstrom^2'), symmetry=1, barrier=(2.83209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0165365,'amu*angstrom^2'), symmetry=1, barrier=(2.83209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0165365,'amu*angstrom^2'), symmetry=1, barrier=(2.83209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0165365,'amu*angstrom^2'), symmetry=1, barrier=(2.83209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0165365,'amu*angstrom^2'), symmetry=1, barrier=(2.83209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02465,0.0713679,-8.81422e-05,7.41875e-08,-2.62932e-11,80253.4,35.3692], Tmin=(100,'K'), Tmax=(816.189,'K')), NASAPolynomial(coeffs=[4.7125,0.043652,-1.94845e-05,3.63305e-09,-2.48681e-13,79972.6,20.2941], Tmin=(816.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(RCCJCC) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])C[CH][CH]C(4457)',
    structure = SMILES('[CH2][C]([O])C[CH][CH]C'),
    E0 = (643.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,360,370,350,3000,3100,440,815,1455,1000,318.743,1982.27,2173.91,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0932087,'amu*angstrom^2'), symmetry=1, barrier=(5.76288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0932087,'amu*angstrom^2'), symmetry=1, barrier=(5.76288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0932087,'amu*angstrom^2'), symmetry=1, barrier=(5.76288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0932087,'amu*angstrom^2'), symmetry=1, barrier=(5.76288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0932087,'amu*angstrom^2'), symmetry=1, barrier=(5.76288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6413,0.0857536,-0.000135476,1.28197e-07,-4.67632e-11,77462.3,33.671], Tmin=(100,'K'), Tmax=(863.944,'K')), NASAPolynomial(coeffs=[2.97831,0.0476613,-2.19881e-05,4.08538e-09,-2.75807e-13,78076.3,28.6285], Tmin=(863.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(C2CsJOH) + radical(CJCO)"""),
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
    E0 = (671.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (837.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (829.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1145.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1321.33,'kJ/mol'),
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
    E0 = (1137.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1137.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1120.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1126.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (677.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (680.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (680.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (694.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (694.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (735.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (735.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (735.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (681.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (686.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (783.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (717.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (679.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (720.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1023.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1083.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1060.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (814.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (785.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (813.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (787.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (747.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (789.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (851.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (824.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (747.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (745.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['vinoxy(1351)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C([CH2])[O](1366)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2][CH][CH]C[CH][CH2](653)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C[CH][O](1850)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][CH][C]CC([CH2])[O](9315)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C][CH]CC([CH2])[O](9316)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([O])C[CH][CH][CH2](9018)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH][CH][CH]CC([CH2])[O](9317)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]CC1CO1(9318)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH]C1CC([CH2])O1(5565)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH]C1CC([O])C1(5529)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]CC(=C)O(4326)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]CC(C)=O(4331)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2]C=C[CH]C([CH2])O(5620)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2]C=C[CH]C(C)[O](5621)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH]CCC(=C)[O](4322)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2]C=C[CH]C([CH2])[O](5434)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][CH][CH]CC(=C)[O](4319)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(T)(20)', '[CH2]C=CC[CH][O](1848)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['vinoxy(1351)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4679.9,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C([CH2])[O](9139)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][CH][CH]C[C]([CH2])[O](9319)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]C[CH]C([CH2])[O](5614)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]C[C]([CH2])O(9320)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]C[C](C)[O](9321)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH]CC[C]([CH2])[O](5615)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH][CH]C([CH2])O(9322)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH][CH]C(C)[O](9323)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C[CH][CH]C([CH2])[O](9324)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2]C[CH]C[C]([CH2])[O](5616)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_1;Y_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2]C([O])[CH][CH][CH]C(4312)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.32767e+06,'s^-1'), n=1.53625, Ea=(75.6258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][C]([O])C[CH][CH]C(4457)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2099',
    isomers = [
        '[CH2][CH][CH]CC([CH2])[O](4458)',
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
    label = 'PDepNetwork #2099',
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


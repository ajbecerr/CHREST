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
    label = 'allene(458)',
    structure = SMILES('C=C=C'),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.0070462,2.78306e-05,-3.99445e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.705,'K')), NASAPolynomial(coeffs=[6.79956,0.00959979,-3.02068e-06,5.37827e-10,-3.92606e-14,19772.3,-12.7582], Tmin=(949.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH]C([CH2])([CH2])C=C(16643)',
    structure = SMILES('[CH2][CH]C([CH2])([CH2])C=C'),
    E0 = (708.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,180,1177.47],'cm^-1')),
        HinderedRotor(inertia=(0.00306004,'amu*angstrom^2'), symmetry=1, barrier=(2.99108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13284,'amu*angstrom^2'), symmetry=1, barrier=(3.05426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132366,'amu*angstrom^2'), symmetry=1, barrier=(3.04337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130167,'amu*angstrom^2'), symmetry=1, barrier=(2.99281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931728,'amu*angstrom^2'), symmetry=1, barrier=(21.4223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.956057,0.0655147,-4.52334e-05,1.59138e-08,-2.32626e-12,85307.6,32.502], Tmin=(100,'K'), Tmax=(1538.41,'K')), NASAPolynomial(coeffs=[12.9265,0.0343906,-1.48863e-05,2.76297e-09,-1.89169e-13,81624.5,-30.4045], Tmin=(1538.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Neopentyl) + radical(Neopentyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])=C(15888)',
    structure = SMILES('[CH2][CH]CC=C([CH2])[CH2]'),
    E0 = (587.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,291.554,1474.19],'cm^-1')),
        HinderedRotor(inertia=(0.00518582,'amu*angstrom^2'), symmetry=1, barrier=(7.99749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365072,'amu*angstrom^2'), symmetry=1, barrier=(22.0206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00518582,'amu*angstrom^2'), symmetry=1, barrier=(7.99749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0354252,'amu*angstrom^2'), symmetry=1, barrier=(54.6308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98317,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756853,0.0664382,-4.59339e-05,1.67474e-08,-2.54131e-12,70757.4,31.612], Tmin=(100,'K'), Tmax=(1505.17,'K')), NASAPolynomial(coeffs=[12.967,0.0339898,-1.35972e-05,2.42502e-09,-1.62456e-13,67081.7,-32.288], Tmin=(1505.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C][CH2](15952)',
    structure = SMILES('[CH2][C][CH2]'),
    E0 = (738.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,618.688],'cm^-1')),
        HinderedRotor(inertia=(0.00838111,'amu*angstrom^2'), symmetry=1, barrier=(2.27602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989004,'amu*angstrom^2'), symmetry=1, barrier=(2.27392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05819,0.0178835,-7.35728e-06,-7.95768e-10,1.01791e-12,88882.1,14.263], Tmin=(100,'K'), Tmax=(1175.48,'K')), NASAPolynomial(coeffs=[6.73198,0.0102413,-3.80615e-06,7.07005e-10,-4.9648e-14,87682.7,-5.48288], Tmin=(1175.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH][CH]C([CH2])=C(15736)',
    structure = SMILES('[CH2][CH]C=C([CH2])[CH2]'),
    E0 = (557.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,933.961],'cm^-1')),
        HinderedRotor(inertia=(0.0997163,'amu*angstrom^2'), symmetry=1, barrier=(110.952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0210093,'amu*angstrom^2'), symmetry=1, barrier=(23.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207712,'amu*angstrom^2'), symmetry=1, barrier=(23.2561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.84831,'amu*angstrom^2'), symmetry=1, barrier=(111.472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38143,0.046517,-5.02403e-06,-2.60325e-08,1.3016e-11,67186.5,23.9745], Tmin=(100,'K'), Tmax=(1001.85,'K')), NASAPolynomial(coeffs=[13.1761,0.0244101,-9.33281e-06,1.72725e-09,-1.22825e-13,63569.3,-39.2074], Tmin=(1001.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P) + radical(Allyl_P) + radical(RCCJ)"""),
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
    label = '[CH]C(C=C)[C]([CH2])[CH2](16644)',
    structure = SMILES('[CH]C(C=C)[C]([CH2])[CH2]'),
    E0 = (942.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,257.455,872.507,1163.34,1524.85,1877.48],'cm^-1')),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75675,0.073724,-8.24397e-05,5.77774e-08,-1.69486e-11,113444,32.9392], Tmin=(100,'K'), Tmax=(865.264,'K')), NASAPolynomial(coeffs=[8.33664,0.0365713,-1.43717e-05,2.51177e-09,-1.65715e-13,112211,-2.07554], Tmin=(865.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(942.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]([CH2])C([CH2])C=C(16645)',
    structure = SMILES('[CH][C]([CH2])C([CH2])C=C'),
    E0 = (942.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,257.455,872.507,1163.34,1524.85,1877.48],'cm^-1')),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086489,'amu*angstrom^2'), symmetry=1, barrier=(3.29519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75675,0.073724,-8.24397e-05,5.77774e-08,-1.69486e-11,113444,33.6324], Tmin=(100,'K'), Tmax=(865.264,'K')), NASAPolynomial(coeffs=[8.33664,0.0365713,-1.43717e-05,2.51177e-09,-1.65715e-13,112211,-1.3824], Tmin=(865.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(942.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1([CH2])CC1C=C(15854)',
    structure = SMILES('[CH2]C1([CH2])CC1C=C'),
    E0 = (455.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10232,0.0448963,3.25721e-05,-7.58248e-08,3.24645e-11,54924.1,26.9751], Tmin=(100,'K'), Tmax=(970.103,'K')), NASAPolynomial(coeffs=[16.932,0.0256963,-8.97509e-06,1.68021e-09,-1.24393e-13,49685,-60.0865], Tmin=(970.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)C=C(5210)',
    structure = SMILES('[CH2]C([CH2])C(=C)C=C'),
    E0 = (412.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,255.308,255.309],'cm^-1')),
        HinderedRotor(inertia=(0.00258603,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00258629,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453229,'amu*angstrom^2'), symmetry=1, barrier=(20.9652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60557,'amu*angstrom^2'), symmetry=1, barrier=(74.267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.573762,0.0643475,-2.34154e-05,-2.22559e-08,1.59976e-11,49697.5,28.677], Tmin=(100,'K'), Tmax=(912.46,'K')), NASAPolynomial(coeffs=[16.2624,0.0249727,-7.0187e-06,1.07672e-09,-7.02751e-14,45610.5,-52.2813], Tmin=(912.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)[C](C)C=C(15859)',
    structure = SMILES('[CH2]C([CH2])=C(C)C=C'),
    E0 = (217.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,854.237],'cm^-1')),
        HinderedRotor(inertia=(0.154235,'amu*angstrom^2'), symmetry=1, barrier=(77.8778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914214,'amu*angstrom^2'), symmetry=1, barrier=(21.0196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.38411,'amu*angstrom^2'), symmetry=1, barrier=(77.8075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383489,'amu*angstrom^2'), symmetry=1, barrier=(77.9394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902769,0.054775,1.97074e-06,-4.17001e-08,2.02939e-11,26266.1,25.3908], Tmin=(100,'K'), Tmax=(959.926,'K')), NASAPolynomial(coeffs=[14.4133,0.0306145,-1.04948e-05,1.83437e-09,-1.27389e-13,22191.6,-46.9497], Tmin=(959.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])[C]([CH2])[CH2](16646)',
    structure = SMILES('[CH2]C[C]([CH2])[C]([CH2])[CH2]'),
    E0 = (967.652,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,363.333,366.667,370,300,400,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00568719,'amu*angstrom^2'), symmetry=1, barrier=(35.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00568719,'amu*angstrom^2'), symmetry=1, barrier=(35.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00568719,'amu*angstrom^2'), symmetry=1, barrier=(35.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00568719,'amu*angstrom^2'), symmetry=1, barrier=(35.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00568719,'amu*angstrom^2'), symmetry=1, barrier=(35.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00568719,'amu*angstrom^2'), symmetry=1, barrier=(35.3243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615745,0.0889246,-0.000145213,1.38624e-07,-4.9362e-11,116490,36.6334], Tmin=(100,'K'), Tmax=(911.336,'K')), NASAPolynomial(coeffs=[1.09123,0.0501775,-2.10972e-05,3.68989e-09,-2.37773e-13,117926,42.736], Tmin=(911.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(967.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Tertalkyl) + radical(RCCJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]CC1([CH2])[CH2](16617)',
    structure = SMILES('[CH2]C1[CH]CC1([CH2])[CH2]'),
    E0 = (715.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26046,0.0462577,1.6916e-05,-5.09962e-08,2.17681e-11,86142.5,30.0424], Tmin=(100,'K'), Tmax=(991.059,'K')), NASAPolynomial(coeffs=[12.9143,0.0323236,-1.20948e-05,2.22052e-09,-1.57492e-13,82206.9,-34.2774], Tmin=(991.059,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(715.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Neopentyl) + radical(Neopentyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C1[CH]CC1(6368)',
    structure = SMILES('[CH2][C]([CH2])C1[CH]CC1'),
    E0 = (707.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,224.77,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,1141.38,2299.12,2299.12],'cm^-1')),
        HinderedRotor(inertia=(0.0881622,'amu*angstrom^2'), symmetry=1, barrier=(2.51912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881622,'amu*angstrom^2'), symmetry=1, barrier=(2.51912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881622,'amu*angstrom^2'), symmetry=1, barrier=(2.51912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34728,0.0517254,-2.01156e-05,2.35035e-10,1.29084e-12,85238.6,31.2295], Tmin=(100,'K'), Tmax=(1279.59,'K')), NASAPolynomial(coeffs=[8.25764,0.0385785,-1.46155e-05,2.53326e-09,-1.67051e-13,82777.9,-6.51744], Tmin=(1279.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(707.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1CC[CH]C1[CH2](15842)',
    structure = SMILES('[CH2][C]1CC[CH]C1[CH2]'),
    E0 = (632.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,301.693,807.677,807.677,807.677,807.677,807.677,807.677,807.677,807.677,807.677,1618.05,1618.05,1618.05,1618.05,1618.05,1618.05,1618.05,1618.05,1618.05],'cm^-1')),
        HinderedRotor(inertia=(0.153016,'amu*angstrom^2'), symmetry=1, barrier=(3.56589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153016,'amu*angstrom^2'), symmetry=1, barrier=(3.56589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86963,0.0371193,2.39651e-05,-4.78665e-08,1.89722e-11,76195,30.4382], Tmin=(100,'K'), Tmax=(976.964,'K')), NASAPolynomial(coeffs=[7.278,0.0383566,-1.38326e-05,2.4226e-09,-1.65061e-13,74022.5,-1.23852], Tmin=(976.964,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C([CH2])C1([CH2])[CH2](16647)',
    structure = SMILES('[CH2]C1C([CH2])C1([CH2])[CH2]'),
    E0 = (727.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896578,0.0531629,8.96942e-06,-5.48347e-08,2.69563e-11,87675.1,29.1027], Tmin=(100,'K'), Tmax=(929.481,'K')), NASAPolynomial(coeffs=[16.5027,0.0245002,-6.90276e-06,1.1107e-09,-7.66459e-14,83111,-53.9923], Tmin=(929.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Neopentyl) + radical(Neopentyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C1CC1[CH2](16648)',
    structure = SMILES('[CH2][C]([CH2])C1CC1[CH2]'),
    E0 = (720.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15461,0.0567438,-2.21172e-05,-1.03243e-08,8.90642e-12,86763.6,30.3602], Tmin=(100,'K'), Tmax=(899.279,'K')), NASAPolynomial(coeffs=[10.0472,0.0336395,-1.10177e-05,1.7884e-09,-1.15766e-13,84499.1,-15.2956], Tmin=(899.279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(720.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1CC([CH2])C1[CH2](16649)',
    structure = SMILES('[CH2][C]1CC([CH2])C1[CH2]'),
    E0 = (716.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41866,0.0480255,2.98591e-06,-3.44097e-08,1.6669e-11,86225.2,30.1182], Tmin=(100,'K'), Tmax=(926.25,'K')), NASAPolynomial(coeffs=[9.65301,0.0347015,-1.14465e-05,1.896e-09,-1.25473e-13,83745.9,-14.126], Tmin=(926.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(716.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])=C(15857)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])C=C'),
    E0 = (368.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,591.3],'cm^-1')),
        HinderedRotor(inertia=(0.219401,'amu*angstrom^2'), symmetry=1, barrier=(54.2451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219624,'amu*angstrom^2'), symmetry=1, barrier=(54.2744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222477,'amu*angstrom^2'), symmetry=1, barrier=(54.2509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221175,'amu*angstrom^2'), symmetry=1, barrier=(54.2562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918381,0.052082,1.12934e-05,-5.6203e-08,2.68704e-11,44489,25.6158], Tmin=(100,'K'), Tmax=(942.248,'K')), NASAPolynomial(coeffs=[16.5778,0.024779,-7.6046e-06,1.29112e-09,-9.1042e-14,39799,-58.2284], Tmin=(942.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])=C(15735)',
    structure = SMILES('[CH2]C([CH2])=CC=C'),
    E0 = (254.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,517.551],'cm^-1')),
        HinderedRotor(inertia=(0.282927,'amu*angstrom^2'), symmetry=1, barrier=(53.9901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281042,'amu*angstrom^2'), symmetry=1, barrier=(53.9846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280507,'amu*angstrom^2'), symmetry=1, barrier=(53.9689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74832,0.0347423,3.26759e-05,-6.90546e-08,2.98664e-11,30758.3,21.7557], Tmin=(100,'K'), Tmax=(938.86,'K')), NASAPolynomial(coeffs=[13.4798,0.0220376,-6.58312e-06,1.11304e-09,-7.90227e-14,26912.5,-42.8515], Tmin=(938.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
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
    label = '[CH2][C]([CH2])C=C(15109)',
    structure = SMILES('[CH2]C=C([CH2])[CH2]'),
    E0 = (385.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.110043,'amu*angstrom^2'), symmetry=1, barrier=(31.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110116,'amu*angstrom^2'), symmetry=1, barrier=(31.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242242,'amu*angstrom^2'), symmetry=1, barrier=(68.3077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82336,0.036309,8.59652e-06,-3.90391e-08,1.80728e-11,46458.6,18.3828], Tmin=(100,'K'), Tmax=(965.965,'K')), NASAPolynomial(coeffs=[13.046,0.0178172,-6.13813e-06,1.11719e-09,-8.09237e-14,42985,-42.1286], Tmin=(965.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][C]([CH2])[CH2](16650)',
    structure = SMILES('[CH2][CH][C]([CH2])[CH2]'),
    E0 = (821.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1649.2],'cm^-1')),
        HinderedRotor(inertia=(0.00214901,'amu*angstrom^2'), symmetry=1, barrier=(4.16088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181089,'amu*angstrom^2'), symmetry=1, barrier=(4.16359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00215337,'amu*angstrom^2'), symmetry=1, barrier=(4.15601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214011,'amu*angstrom^2'), symmetry=1, barrier=(4.13915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07932,0.0503958,-7.54465e-05,7.29853e-08,-2.67542e-11,98876.7,26.8107], Tmin=(100,'K'), Tmax=(899.507,'K')), NASAPolynomial(coeffs=[1.11254,0.0342709,-1.4498e-05,2.57078e-09,-1.68211e-13,99876.9,35.9654], Tmin=(899.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(821.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH2])[C]([CH2])C=C(16651)',
    structure = SMILES('[CH2][CH]C([CH2])=C([CH2])[CH2]'),
    E0 = (670.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,1261.77],'cm^-1')),
        HinderedRotor(inertia=(0.0738167,'amu*angstrom^2'), symmetry=1, barrier=(83.3766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0230441,'amu*angstrom^2'), symmetry=1, barrier=(26.0332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0230379,'amu*angstrom^2'), symmetry=1, barrier=(26.0334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.6261,'amu*angstrom^2'), symmetry=1, barrier=(83.3712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698603,'amu*angstrom^2'), symmetry=1, barrier=(26.0313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.54049,0.0639845,-2.68411e-05,-1.26464e-08,9.80904e-12,80741.5,28.0302], Tmin=(100,'K'), Tmax=(1006.63,'K')), NASAPolynomial(coeffs=[16.3395,0.0270415,-1.02914e-05,1.89055e-09,-1.33624e-13,76251.7,-54.7974], Tmin=(1006.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_S) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH2])C([CH2])[C]=C(16652)',
    structure = SMILES('[CH2][C]([CH2])C([CH2])[C]=C'),
    E0 = (936.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,1685,370,180,2325],'cm^-1')),
        HinderedRotor(inertia=(0.000888784,'amu*angstrom^2'), symmetry=1, barrier=(3.41015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148531,'amu*angstrom^2'), symmetry=1, barrier=(3.41502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148234,'amu*angstrom^2'), symmetry=1, barrier=(3.4082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0175587,'amu*angstrom^2'), symmetry=1, barrier=(67.3623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.92763,'amu*angstrom^2'), symmetry=1, barrier=(67.3119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729135,0.0777491,-0.000103106,8.37152e-08,-2.68029e-11,112805,34.2873], Tmin=(100,'K'), Tmax=(928.593,'K')), NASAPolynomial(coeffs=[6.35596,0.0382338,-1.45966e-05,2.45466e-09,-1.55709e-13,112419,11.1046], Tmin=(928.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(936.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])[C]([CH2])[CH2](14010)',
    structure = SMILES('[CH]=CC([CH2])[C]([CH2])[CH2]'),
    E0 = (946.237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3120,650,792.5,1650,1303.84],'cm^-1')),
        HinderedRotor(inertia=(0.146359,'amu*angstrom^2'), symmetry=1, barrier=(3.36508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146495,'amu*angstrom^2'), symmetry=1, barrier=(3.36821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146404,'amu*angstrom^2'), symmetry=1, barrier=(3.36611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146355,'amu*angstrom^2'), symmetry=1, barrier=(3.36499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00278618,'amu*angstrom^2'), symmetry=1, barrier=(3.36311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1464,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.67329,0.0768515,-9.56014e-05,7.26193e-08,-2.20728e-11,113922,34.417], Tmin=(100,'K'), Tmax=(938.284,'K')), NASAPolynomial(coeffs=[7.84699,0.035796,-1.32241e-05,2.19248e-09,-1.38316e-13,113037,2.72173], Tmin=(938.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(946.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH2][CH][C](C)C([CH2])=C(16035)',
    structure = SMILES('[CH2][CH]C(C)=C([CH2])[CH2]'),
    E0 = (518.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0142132,'amu*angstrom^2'), symmetry=1, barrier=(6.46654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164808,'amu*angstrom^2'), symmetry=1, barrier=(22.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.971863,'amu*angstrom^2'), symmetry=1, barrier=(22.345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0647866,'amu*angstrom^2'), symmetry=1, barrier=(87.856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0268779,'amu*angstrom^2'), symmetry=1, barrier=(36.4322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483806,0.0671015,-3.73251e-05,2.82769e-09,3.08234e-12,62520.5,27.9561], Tmin=(100,'K'), Tmax=(1086.6,'K')), NASAPolynomial(coeffs=[14.694,0.0320295,-1.2707e-05,2.32411e-09,-1.61022e-13,58414.6,-46.4629], Tmin=(1086.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_S) + radical(Allyl_P) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])[C]=C(5214)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[C]=C'),
    E0 = (751.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,1685,370,1131.4,1132.04],'cm^-1')),
        HinderedRotor(inertia=(0.112624,'amu*angstrom^2'), symmetry=1, barrier=(2.58945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112832,'amu*angstrom^2'), symmetry=1, barrier=(2.59423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112607,'amu*angstrom^2'), symmetry=1, barrier=(2.58905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112699,'amu*angstrom^2'), symmetry=1, barrier=(2.59117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.77263,'amu*angstrom^2'), symmetry=1, barrier=(63.7482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705415,0.0694993,-5.28305e-05,1.64668e-08,8.97313e-13,90513.3,33.2435], Tmin=(100,'K'), Tmax=(853.833,'K')), NASAPolynomial(coeffs=[11.4984,0.0322036,-1.06167e-05,1.70421e-09,-1.08476e-13,88186.7,-19.9526], Tmin=(853.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]C(C)C([CH2])=C(16026)',
    structure = SMILES('[CH2][C]C(C)C([CH2])=C'),
    E0 = (698.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.341483,0.0703217,-4.40636e-05,7.04796e-09,2.3356e-12,84094.7,29.4912], Tmin=(100,'K'), Tmax=(1060.23,'K')), NASAPolynomial(coeffs=[15.6371,0.0299988,-1.161e-05,2.10629e-09,-1.4578e-13,79874.3,-49.8037], Tmin=(1060.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(698.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Allyl_P) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C](C)[C]([CH2])C=C(16653)',
    structure = SMILES('[CH2][CH]C([CH2])=C([CH2])C'),
    E0 = (518.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483806,0.0671015,-3.73251e-05,2.82769e-09,3.08234e-12,62520.5,27.9561], Tmin=(100,'K'), Tmax=(1086.6,'K')), NASAPolynomial(coeffs=[14.694,0.0320295,-1.2707e-05,2.32411e-09,-1.61022e-13,58414.6,-46.4629], Tmin=(1086.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_S) + radical(Allyl_P) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C)C([CH2])[C]=C(15848)',
    structure = SMILES('[CH2][C](C)C([CH2])[C]=C'),
    E0 = (731.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,180,1803.71],'cm^-1')),
        HinderedRotor(inertia=(0.0857287,'amu*angstrom^2'), symmetry=1, barrier=(1.97107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0891937,'amu*angstrom^2'), symmetry=1, barrier=(2.05074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859554,'amu*angstrom^2'), symmetry=1, barrier=(1.97628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0818946,'amu*angstrom^2'), symmetry=1, barrier=(1.88292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0808003,'amu*angstrom^2'), symmetry=1, barrier=(1.85776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.775228,0.0752946,-8.7894e-05,6.77765e-08,-2.18071e-11,88139.2,33.0348], Tmin=(100,'K'), Tmax=(871.73,'K')), NASAPolynomial(coeffs=[6.24434,0.0421593,-1.70432e-05,3.01232e-09,-1.99414e-13,87491.2,9.15247], Tmin=(871.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(731.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C([CH2])[CH2](5215)',
    structure = SMILES('[CH]=CC([CH2])C([CH2])[CH2]'),
    E0 = (760.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3120,650,792.5,1650,3054.54],'cm^-1')),
        HinderedRotor(inertia=(0.0118335,'amu*angstrom^2'), symmetry=1, barrier=(78.3593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00261145,'amu*angstrom^2'), symmetry=1, barrier=(17.3057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40781,'amu*angstrom^2'), symmetry=1, barrier=(78.3523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40772,'amu*angstrom^2'), symmetry=1, barrier=(78.3503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40913,'amu*angstrom^2'), symmetry=1, barrier=(78.3827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3397.52,'J/mol'), sigma=(6.23225,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.69 K, Pc=31.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.182796,0.0742221,-6.55739e-05,3.25935e-08,-6.47757e-12,91650.5,35.0403], Tmin=(100,'K'), Tmax=(1311.98,'K')), NASAPolynomial(coeffs=[14.1921,0.0277129,-8.0581e-06,1.16155e-09,-6.78052e-14,88301.3,-35.1053], Tmin=(1311.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(C)[C]([CH2])[CH2](14013)',
    structure = SMILES('[CH]=CC(C)[C]([CH2])[CH2]'),
    E0 = (741.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.219409,'amu*angstrom^2'), symmetry=1, barrier=(5.04464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219106,'amu*angstrom^2'), symmetry=1, barrier=(5.03768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87425,'amu*angstrom^2'), symmetry=1, barrier=(66.0846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186199,'amu*angstrom^2'), symmetry=1, barrier=(5.03861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186036,'amu*angstrom^2'), symmetry=1, barrier=(5.03845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.750146,0.0740073,-7.88943e-05,5.45373e-08,-1.6065e-11,89254.9,32.3627], Tmin=(100,'K'), Tmax=(852.844,'K')), NASAPolynomial(coeffs=[7.65689,0.0398605,-1.57533e-05,2.77006e-09,-1.83699e-13,88140.6,0.514806], Tmin=(852.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])[C]([CH2])C(13826)',
    structure = SMILES('[CH]=CC([CH2])[C]([CH2])C'),
    E0 = (741.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.219409,'amu*angstrom^2'), symmetry=1, barrier=(5.04464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219106,'amu*angstrom^2'), symmetry=1, barrier=(5.03768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87425,'amu*angstrom^2'), symmetry=1, barrier=(66.0846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186199,'amu*angstrom^2'), symmetry=1, barrier=(5.03861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186036,'amu*angstrom^2'), symmetry=1, barrier=(5.03845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.750146,0.0740073,-7.88943e-05,5.45373e-08,-1.6065e-11,89254.9,33.0559], Tmin=(100,'K'), Tmax=(852.844,'K')), NASAPolynomial(coeffs=[7.65689,0.0398605,-1.57533e-05,2.77006e-09,-1.83699e-13,88140.6,1.20795], Tmin=(852.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
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
    E0 = (699.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (871.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (856.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1068.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (995.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1154.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1154.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (704.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (722.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (722.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (990.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (838.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (881.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (819.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (754.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (736.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (792.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (699.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (699.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (851.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (699.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1108.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (881.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1148.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1158.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (857.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (840.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (893.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (849.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (856.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (835.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (805.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (785.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (774.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['allene(458)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])([CH2])C=C(16643)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][CH]C[CH]C([CH2])=C(15888)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][C][CH2](15952)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C([CH2])=C(15736)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C(C=C)[C]([CH2])[CH2](16644)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH][C]([CH2])C([CH2])C=C(16645)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2]C1([CH2])CC1C=C(15854)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2]C([CH2])C(=C)C=C(5210)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2]C(=C)[C](C)C=C(15859)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C[C]([CH2])[C]([CH2])[CH2](16646)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2]C1[CH]CC1([CH2])[CH2](16617)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.61918e+08,'s^-1'), n=0.930343, Ea=(139.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][C]([CH2])C1[CH]CC1(6368)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.18363e+11,'s^-1'), n=0.209288, Ea=(182.141,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 125 used for R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][C]1CC[CH]C1[CH2](15842)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4779e+10,'s^-1'), n=0.565913, Ea=(120.257,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2]C1C([CH2])C1([CH2])[CH2](16647)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.65211e+08,'s^-1'), n=0.77, Ea=(55.5948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][C]([CH2])C1CC1[CH2](16648)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][C]1CC([CH2])C1[CH2](16649)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.24689e+08,'s^-1'), n=0.9705, Ea=(92.9036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][C](C=C)C([CH2])=C(15857)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(118.488,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 114.7 to 118.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(T)(20)', '[CH2]C=CC([CH2])=C(15735)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(61.777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 57.4 to 61.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C][CH2](15952)', 'butadiene13(1350)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.10641,'m^3/(mol*s)'), n=2.065, Ea=(15.9938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;YJ] for rate rule [Cds-CdH_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C2H3(60)', '[CH2][C]([CH2])C=C(15109)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0131003,'m^3/(mol*s)'), n=2.40999, Ea=(26.6398,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 22.1 to 26.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C2H3(60)', '[CH2][CH][C]([CH2])[CH2](16650)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.53609e+09,'m^3/(mol*s)'), n=-0.946459, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.000724761308015, var=1.22133425516, Tref=1000.0, N=6, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_Ext-2R-R_Ext-1C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_Ext-2R-R_Ext-1C-R_6R!H->C
    Total Standard Deviation in ln(k): 2.2173337826
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_Ext-2R-R_Ext-1C-R_6R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C]([CH2])[C]([CH2])C=C(16651)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][C]([CH2])C([CH2])[C]=C(16652)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]=CC([CH2])[C]([CH2])[CH2](14010)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][C](C=C)C([CH2])[CH2](5212)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][CH][C](C)C([CH2])=C(16035)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])C([CH2])[C]=C(5214)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]C(C)C([CH2])=C(16026)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    products = ['[CH2][C](C)[C]([CH2])C=C(16653)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(19639.3,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C](C)C([CH2])[C]=C(15848)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CC([CH2])C([CH2])[CH2](5215)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC(C)[C]([CH2])[CH2](14013)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=CC([CH2])[C]([CH2])C(13826)'],
    products = ['[CH2][C]([CH2])C([CH2])C=C(5213)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3655',
    isomers = [
        '[CH2][C]([CH2])C([CH2])C=C(5213)',
    ],
    reactants = [
        ('allene(458)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3655',
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


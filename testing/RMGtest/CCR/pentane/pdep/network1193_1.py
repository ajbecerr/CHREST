species(
    label = '[CH2]C([CH]C)C(=C)[O](2503)',
    structure = SMILES('[CH2]C([CH]C)C(=C)[O]'),
    E0 = (250.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,1962.86],'cm^-1')),
        HinderedRotor(inertia=(0.0102752,'amu*angstrom^2'), symmetry=1, barrier=(2.7886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0102751,'amu*angstrom^2'), symmetry=1, barrier=(2.78855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373544,'amu*angstrom^2'), symmetry=1, barrier=(10.1375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373544,'amu*angstrom^2'), symmetry=1, barrier=(10.1375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586634,0.0664197,-5.34853e-05,2.33575e-08,-4.13471e-12,30296.4,32.4745], Tmin=(100,'K'), Tmax=(1352.09,'K')), NASAPolynomial(coeffs=[14.0751,0.0265156,-9.2158e-06,1.52979e-09,-9.87714e-14,26648.9,-36.6682], Tmin=(1352.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'ketene(1375)',
    structure = SMILES('C=C=O'),
    E0 = (-60.9858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48302,0.0083299,5.93417e-06,-1.33939e-08,5.73133e-12,-7313.53,6.51898], Tmin=(100,'K'), Tmax=(954.863,'K')), NASAPolynomial(coeffs=[5.88246,0.00580173,-1.91266e-06,3.35917e-10,-2.37398e-14,-8114.73,-6.74202], Tmin=(954.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH]C(C)C(=C)[O](2743)',
    structure = SMILES('[CH2][CH]C(C)C(=C)[O]'),
    E0 = (250.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,997.736,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0268206,'amu*angstrom^2'), symmetry=1, barrier=(18.9471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00739858,'amu*angstrom^2'), symmetry=1, barrier=(5.22675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227331,'amu*angstrom^2'), symmetry=1, barrier=(5.22678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.824074,'amu*angstrom^2'), symmetry=1, barrier=(18.9471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.58,'J/mol'), sigma=(6.61681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.14 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645949,0.0649023,-4.90165e-05,1.93488e-08,-3.09336e-12,30314.1,32.2727], Tmin=(100,'K'), Tmax=(1479.08,'K')), NASAPolynomial(coeffs=[14.9066,0.0263364,-9.90575e-06,1.7206e-09,-1.13823e-13,26095.5,-42.1091], Tmin=(1479.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C[CH][CH]C(2560)',
    structure = SMILES('C=C([O])C[CH][CH]C'),
    E0 = (248.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,222.408,1946.36,1946.36,1946.36],'cm^-1')),
        HinderedRotor(inertia=(0.20821,'amu*angstrom^2'), symmetry=1, barrier=(7.30844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00340802,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208215,'amu*angstrom^2'), symmetry=1, barrier=(7.30845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00340791,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02514,0.0519567,-1.23522e-06,-9.02755e-08,8.70563e-11,29898,28.692], Tmin=(100,'K'), Tmax=(467.454,'K')), NASAPolynomial(coeffs=[4.96563,0.0420012,-1.80839e-05,3.34267e-09,-2.2893e-13,29457,14.9649], Tmin=(467.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC)"""),
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
    label = '[CH2][CH]C(=C)[O](2379)',
    structure = SMILES('[CH2]C=C([CH2])[O]'),
    E0 = (203.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,582.136],'cm^-1')),
        HinderedRotor(inertia=(0.0360605,'amu*angstrom^2'), symmetry=1, barrier=(8.66353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0960342,'amu*angstrom^2'), symmetry=1, barrier=(22.9656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97595,0.036843,-1.21657e-05,-1.70174e-08,1.14342e-11,24586.7,18.7658], Tmin=(100,'K'), Tmax=(913.194,'K')), NASAPolynomial(coeffs=[12.671,0.0104313,-2.34913e-06,3.21192e-10,-2.11198e-14,21781.3,-36.526], Tmin=(913.194,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
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
    label = 'C=C([O])[CH][CH]C(2452)',
    structure = SMILES('[CH2]C([O])=C[CH]C'),
    E0 = (170.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,180.08,607.165],'cm^-1')),
        HinderedRotor(inertia=(0.0075904,'amu*angstrom^2'), symmetry=1, barrier=(1.99596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0878727,'amu*angstrom^2'), symmetry=1, barrier=(2.02037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0699717,'amu*angstrom^2'), symmetry=1, barrier=(18.398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49458,0.0463747,-1.57054e-05,-1.54896e-08,1.03317e-11,20632,22.0175], Tmin=(100,'K'), Tmax=(952.985,'K')), NASAPolynomial(coeffs=[12.8894,0.0195111,-6.42011e-06,1.09879e-09,-7.56641e-14,17508.2,-37.402], Tmin=(952.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=C(O)CJ)"""),
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
    label = '[CH]C([CH2])C(=C)[O](2593)',
    structure = SMILES('[CH]C([CH2])C(=C)[O]'),
    E0 = (528.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,330.312,330.546,330.591,1084.67,1084.73],'cm^-1')),
        HinderedRotor(inertia=(0.00154243,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000143262,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110272,'amu*angstrom^2'), symmetry=1, barrier=(8.58286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803084,0.0581349,-5.56878e-05,2.73487e-08,-5.12898e-12,63661.7,27.8988], Tmin=(100,'K'), Tmax=(1475.95,'K')), NASAPolynomial(coeffs=[14.948,0.0129381,-2.78021e-06,3.0083e-10,-1.39499e-14,60233.7,-43.3167], Tmin=(1475.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C([C]=C)[CH]C(4287)',
    structure = SMILES('[CH2]C([C]=C)[CH]C'),
    E0 = (565.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180,823.069],'cm^-1')),
        HinderedRotor(inertia=(0.13717,'amu*angstrom^2'), symmetry=1, barrier=(3.1538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133446,'amu*angstrom^2'), symmetry=1, barrier=(3.06819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13228,'amu*angstrom^2'), symmetry=1, barrier=(3.04138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0065176,'amu*angstrom^2'), symmetry=1, barrier=(3.16416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57205,0.0518488,-3.31338e-05,1.15546e-08,-1.73452e-12,68091.7,27.8445], Tmin=(100,'K'), Tmax=(1469.6,'K')), NASAPolynomial(coeffs=[8.96437,0.0317281,-1.2597e-05,2.23834e-09,-1.49685e-13,65919,-10.6651], Tmin=(1469.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(Cds_S)"""),
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
    label = '[CH2]C([C]C)C(=C)[O](4288)',
    structure = SMILES('[CH2]C([C]C)C(=C)[O]'),
    E0 = (504.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,266.036,706.615,707.598],'cm^-1')),
        HinderedRotor(inertia=(0.642969,'amu*angstrom^2'), symmetry=1, barrier=(14.7831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.642996,'amu*angstrom^2'), symmetry=1, barrier=(14.7837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0026705,'amu*angstrom^2'), symmetry=1, barrier=(14.7798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221955,'amu*angstrom^2'), symmetry=1, barrier=(78.8095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.198904,0.0724378,-6.79198e-05,3.33775e-08,-6.4168e-12,60822.4,31.2127], Tmin=(100,'K'), Tmax=(1350.48,'K')), NASAPolynomial(coeffs=[16.7738,0.0199406,-5.82941e-06,8.60102e-10,-5.16852e-14,56656,-52.5823], Tmin=(1350.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH]C)C(=C)[O](4289)',
    structure = SMILES('[CH]C([CH]C)C(=C)[O]'),
    E0 = (493.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626373,0.0664894,-5.59305e-05,2.45652e-08,-4.33096e-12,59536.2,31.5544], Tmin=(100,'K'), Tmax=(1359.56,'K')), NASAPolynomial(coeffs=[15.0606,0.0240223,-9.0769e-06,1.59042e-09,-1.06316e-13,55611.4,-42.5162], Tmin=(1359.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C([O])C1CC1C(4290)',
    structure = SMILES('C=C([O])C1CC1C'),
    E0 = (4.74696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30016,0.0389863,4.57645e-05,-9.32769e-08,4.03847e-11,686.801,24.1508], Tmin=(100,'K'), Tmax=(941.233,'K')), NASAPolynomial(coeffs=[18.0417,0.0193564,-5.0527e-06,8.67546e-10,-6.63628e-14,-4746.74,-67.7258], Tmin=(941.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.74696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1C(=C)OC1C(4188)',
    structure = SMILES('[CH2]C1C(=C)OC1C'),
    E0 = (81.6696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09822,0.0410362,5.40254e-05,-1.15495e-07,5.3146e-11,9948.51,20.3217], Tmin=(100,'K'), Tmax=(897.631,'K')), NASAPolynomial(coeffs=[21.3914,0.011624,1.21031e-06,-5.40427e-10,3.86454e-14,3847.11,-89.0827], Tmin=(897.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.6696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Isobutyl)"""),
)

species(
    label = 'C=C1OCC1[CH]C(4273)',
    structure = SMILES('C=C1OCC1[CH]C'),
    E0 = (83.7723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26968,0.0404139,4.34873e-05,-9.35336e-08,4.17635e-11,10192,21.5158], Tmin=(100,'K'), Tmax=(920.736,'K')), NASAPolynomial(coeffs=[17.9998,0.0185296,-3.6161e-06,4.91943e-10,-3.6084e-14,4957.99,-69.5086], Tmin=(920.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.7723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cs_S)"""),
)

species(
    label = 'C=C([O])C(=C)CC(2501)',
    structure = SMILES('C=C([O])C(=C)CC'),
    E0 = (-42.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,285.427,285.702,285.899],'cm^-1')),
        HinderedRotor(inertia=(0.236473,'amu*angstrom^2'), symmetry=1, barrier=(13.688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236176,'amu*angstrom^2'), symmetry=1, barrier=(13.6874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236635,'amu*angstrom^2'), symmetry=1, barrier=(13.6873,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463708,0.0681642,-4.34909e-05,3.06432e-09,5.28816e-12,-4963.98,24.9441], Tmin=(100,'K'), Tmax=(971.204,'K')), NASAPolynomial(coeffs=[16.3402,0.0240564,-8.2358e-06,1.42601e-09,-9.78452e-14,-9051.51,-56.3542], Tmin=(971.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(C)=CC(2740)',
    structure = SMILES('C=C([O])C(C)=CC'),
    E0 = (-55.7831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.574194,'amu*angstrom^2'), symmetry=1, barrier=(13.2019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.573735,'amu*angstrom^2'), symmetry=1, barrier=(13.1913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575063,'amu*angstrom^2'), symmetry=1, barrier=(13.2218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381902,0.0728234,-6.47557e-05,3.07487e-08,-5.86907e-12,-6573.24,24.1332], Tmin=(100,'K'), Tmax=(1262.52,'K')), NASAPolynomial(coeffs=[15.0118,0.0264716,-9.68466e-06,1.66844e-09,-1.10648e-13,-10267.3,-49.8579], Tmin=(1262.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.7831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=CC)C(=C)O(4291)',
    structure = SMILES('[CH2]C(=CC)C(=C)O'),
    E0 = (-42.0887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121784,0.0730736,-4.45292e-05,-4.84468e-09,1.02427e-11,-4911.33,23.845], Tmin=(100,'K'), Tmax=(945.737,'K')), NASAPolynomial(coeffs=[19.9711,0.0188286,-5.61151e-06,9.36085e-10,-6.54554e-14,-9994.33,-77.8328], Tmin=(945.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.0887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC(C)C(=C)[O](2741)',
    structure = SMILES('C=CC(C)C(=C)[O]'),
    E0 = (-21.5737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,326.269,326.427,327.172],'cm^-1')),
        HinderedRotor(inertia=(0.134119,'amu*angstrom^2'), symmetry=1, barrier=(10.1419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133961,'amu*angstrom^2'), symmetry=1, barrier=(10.1419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133922,'amu*angstrom^2'), symmetry=1, barrier=(10.1391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.565508,0.0683287,-5.57221e-05,2.41878e-08,-4.25649e-12,-2465.11,26.2997], Tmin=(100,'K'), Tmax=(1353.73,'K')), NASAPolynomial(coeffs=[14.4619,0.0272675,-1.02241e-05,1.78152e-09,-1.18598e-13,-6227.5,-44.9512], Tmin=(1353.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.5737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(C=C)C(=C)O(3819)',
    structure = SMILES('[CH2]C(C=C)C(=C)O'),
    E0 = (45.7038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0827865,0.0775295,-7.31474e-05,3.58698e-08,-6.86e-12,5654.3,28.9068], Tmin=(100,'K'), Tmax=(1359.35,'K')), NASAPolynomial(coeffs=[18.0989,0.0203041,-5.89138e-06,8.69968e-10,-5.24696e-14,1055.33,-63.1257], Tmin=(1359.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH]C)C([CH2])[O](4292)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])[O]'),
    E0 = (627.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,261.475,954.116,2097.57],'cm^-1')),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903515,0.0667947,-5.41555e-05,2.44318e-08,-4.63745e-12,75579.2,32.26], Tmin=(100,'K'), Tmax=(1227.12,'K')), NASAPolynomial(coeffs=[11.031,0.0337824,-1.38022e-05,2.50879e-09,-1.71098e-13,73093.7,-18.6723], Tmin=(1227.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)[C]1CO1(4293)',
    structure = SMILES('[CH2]C([CH]C)[C]1CO1'),
    E0 = (399.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.795796,0.06709,-6.27193e-05,3.65747e-08,-8.71009e-12,48138.6,30.1576], Tmin=(100,'K'), Tmax=(1136.62,'K')), NASAPolynomial(coeffs=[9.24273,0.0321232,-9.65804e-06,1.39615e-09,-8.03855e-14,46556.9,-10.1867], Tmin=(1136.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[C]([O])CC1C(4294)',
    structure = SMILES('[CH2]C1[C]([O])CC1C'),
    E0 = (387.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20805,0.0473345,1.27246e-05,-5.0522e-08,2.30267e-11,46758.4,25.8006], Tmin=(100,'K'), Tmax=(962.211,'K')), NASAPolynomial(coeffs=[14.2405,0.0273195,-9.32975e-06,1.65666e-09,-1.17118e-13,42668.9,-44.7894], Tmin=(962.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]C1CC[C]1[O](4295)',
    structure = SMILES('C[CH]C1CC[C]1[O]'),
    E0 = (386.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3033,0.0476582,1.57986e-06,-3.00829e-08,1.30679e-11,46572,27.2104], Tmin=(100,'K'), Tmax=(1052.38,'K')), NASAPolynomial(coeffs=[12.0565,0.0324422,-1.33011e-05,2.51001e-09,-1.78093e-13,42888,-31.9666], Tmin=(1052.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C1C(C)C1([CH2])[O](4212)',
    structure = SMILES('[CH2]C1C(C)C1([CH2])[O]'),
    E0 = (408.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596307,0.0626038,-2.25276e-05,-2.16574e-08,1.49012e-11,49229.9,25.7244], Tmin=(100,'K'), Tmax=(944.379,'K')), NASAPolynomial(coeffs=[17.1068,0.022775,-7.07934e-06,1.19017e-09,-8.23406e-14,44769.2,-60.0912], Tmin=(944.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(Isobutyl) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])CC1[CH]C(4296)',
    structure = SMILES('[CH2]C1([O])CC1[CH]C'),
    E0 = (406.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.657555,0.0633488,-3.52462e-05,9.76606e-10,3.93287e-12,49045,27.2548], Tmin=(100,'K'), Tmax=(1062.77,'K')), NASAPolynomial(coeffs=[14.9926,0.0277752,-1.09785e-05,2.02619e-09,-1.41862e-13,44960.1,-47.6593], Tmin=(1062.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(Cs_S) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(=CC)C(=C)[O](4297)',
    structure = SMILES('[CH2]C(=CC)C(=C)[O]'),
    E0 = (95.7161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,220.319,220.561],'cm^-1')),
        HinderedRotor(inertia=(0.501691,'amu*angstrom^2'), symmetry=1, barrier=(17.2542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500434,'amu*angstrom^2'), symmetry=1, barrier=(17.2556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500653,'amu*angstrom^2'), symmetry=1, barrier=(17.2524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462914,0.0694054,-5.31693e-05,1.38006e-08,1.49844e-12,11646.8,24.1211], Tmin=(100,'K'), Tmax=(975.079,'K')), NASAPolynomial(coeffs=[16.3564,0.0219777,-7.54728e-06,1.29949e-09,-8.85418e-14,7702.54,-56.4864], Tmin=(975.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.7161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C(=C)[O](4298)',
    structure = SMILES('[CH2]C(C=C)C(=C)[O]'),
    E0 = (183.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,328.501,328.502,2023.4],'cm^-1')),
        HinderedRotor(inertia=(0.17014,'amu*angstrom^2'), symmetry=1, barrier=(13.0289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170141,'amu*angstrom^2'), symmetry=1, barrier=(13.029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170146,'amu*angstrom^2'), symmetry=1, barrier=(13.0289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648501,0.0693304,-6.62578e-05,3.47965e-08,-7.35884e-12,22195.4,27.7785], Tmin=(100,'K'), Tmax=(1145.95,'K')), NASAPolynomial(coeffs=[13.1664,0.0256346,-9.05998e-06,1.52018e-09,-9.90661e-14,19326.5,-34.3183], Tmin=(1145.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C([O])C=CC(2451)',
    structure = SMILES('C=C([O])C=CC'),
    E0 = (-18.1925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.679055,'amu*angstrom^2'), symmetry=1, barrier=(15.6128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681905,'amu*angstrom^2'), symmetry=1, barrier=(15.6783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.301,0.0519692,-3.14497e-05,5.1886e-10,4.67221e-12,-2084.25,20.2318], Tmin=(100,'K'), Tmax=(967.227,'K')), NASAPolynomial(coeffs=[13.2178,0.0193052,-6.56545e-06,1.13077e-09,-7.7298e-14,-5166.84,-40.8811], Tmin=(967.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.1925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=C)[O](2319)',
    structure = SMILES('C=CC(=C)[O]'),
    E0 = (17.8331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.91425,'amu*angstrom^2'), symmetry=1, barrier=(21.0204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98668,0.0359126,-9.09535e-06,-2.04236e-08,1.25939e-11,2225.07,15.9658], Tmin=(100,'K'), Tmax=(920.926,'K')), NASAPolynomial(coeffs=[13.0167,0.0098341,-2.17549e-06,3.06812e-10,-2.11582e-14,-732.212,-41.3652], Tmin=(920.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=O(1376)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,623.763,623.847],'cm^-1')),
        HinderedRotor(inertia=(0.000434039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44947,0.00889162,4.91274e-06,-1.14995e-08,4.62261e-12,19378.1,9.99239], Tmin=(100,'K'), Tmax=(1030.77,'K')), NASAPolynomial(coeffs=[6.10651,0.00625992,-2.43247e-06,4.78594e-10,-3.54632e-14,18422.4,-4.88573], Tmin=(1030.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJC=O) + radical(CsCJ=O)"""),
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
    label = '[CH2][C]([CH]C)C(=C)[O](4299)',
    structure = SMILES('[CH2]C([O])=C([CH2])[CH]C'),
    E0 = (283.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,601.645,606.348],'cm^-1')),
        HinderedRotor(inertia=(0.0109853,'amu*angstrom^2'), symmetry=1, barrier=(2.82248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0107848,'amu*angstrom^2'), symmetry=1, barrier=(2.87349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828046,'amu*angstrom^2'), symmetry=1, barrier=(19.0384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30746,'amu*angstrom^2'), symmetry=1, barrier=(80.9628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660346,0.0637655,-3.72678e-05,-2.40657e-09,7.23921e-12,34186.7,26.049], Tmin=(100,'K'), Tmax=(959.024,'K')), NASAPolynomial(coeffs=[16.008,0.0222175,-7.42132e-06,1.27208e-09,-8.72857e-14,30209.9,-52.7388], Tmin=(959.024,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(=C)[O](4216)',
    structure = SMILES('[CH2][CH]C([CH2])C(=C)[O]'),
    E0 = (456.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,350.492,1313.27,1313.36],'cm^-1')),
        HinderedRotor(inertia=(0.00133196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00122576,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0519993,'amu*angstrom^2'), symmetry=1, barrier=(4.79708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475557,'amu*angstrom^2'), symmetry=1, barrier=(4.86684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796865,0.065193,-5.74532e-05,2.77245e-08,-5.42681e-12,54971.4,33.5014], Tmin=(100,'K'), Tmax=(1227.73,'K')), NASAPolynomial(coeffs=[12.8885,0.0257971,-9.31944e-06,1.58696e-09,-1.0435e-13,52002.5,-27.3141], Tmin=(1227.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C([O])C([CH2])[CH]C(4300)',
    structure = SMILES('[CH]=C([O])C([CH2])[CH]C'),
    E0 = (497.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,1980.37],'cm^-1')),
        HinderedRotor(inertia=(0.0131234,'amu*angstrom^2'), symmetry=1, barrier=(2.2738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0134473,'amu*angstrom^2'), symmetry=1, barrier=(2.32941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0582614,'amu*angstrom^2'), symmetry=1, barrier=(10.3351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0596179,'amu*angstrom^2'), symmetry=1, barrier=(10.3424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709318,0.0677158,-6.27834e-05,3.18125e-08,-6.50661e-12,60007.5,32.4317], Tmin=(100,'K'), Tmax=(1180.78,'K')), NASAPolynomial(coeffs=[13.1574,0.0255471,-9.21485e-06,1.56793e-09,-1.03116e-13,57067.8,-29.6915], Tmin=(1180.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](CC)C(=C)[O](2502)',
    structure = SMILES('[CH2]C([O])=C([CH2])CC'),
    E0 = (142.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,510.676],'cm^-1')),
        HinderedRotor(inertia=(0.122182,'amu*angstrom^2'), symmetry=1, barrier=(2.80921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0151773,'amu*angstrom^2'), symmetry=1, barrier=(2.80753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0898592,'amu*angstrom^2'), symmetry=1, barrier=(16.6328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0899122,'amu*angstrom^2'), symmetry=1, barrier=(16.6344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451218,0.0691213,-4.66874e-05,6.69044e-09,4.00509e-12,17221.6,27.906], Tmin=(100,'K'), Tmax=(967.351,'K')), NASAPolynomial(coeffs=[15.9798,0.0246775,-8.4227e-06,1.44343e-09,-9.8054e-14,13292.4,-51.2754], Tmin=(967.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])[C](C)[CH]C(2742)',
    structure = SMILES('[CH2]C([O])=C(C)[CH]C'),
    E0 = (131.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,819.675],'cm^-1')),
        HinderedRotor(inertia=(0.0082949,'amu*angstrom^2'), symmetry=1, barrier=(4.05512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60485,'amu*angstrom^2'), symmetry=1, barrier=(13.9067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174309,'amu*angstrom^2'), symmetry=1, barrier=(4.00772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0293364,'amu*angstrom^2'), symmetry=1, barrier=(13.967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.604356,0.0669232,-4.81416e-05,1.39616e-08,-5.16752e-14,15965.6,25.9694], Tmin=(100,'K'), Tmax=(1028.35,'K')), NASAPolynomial(coeffs=[14.1014,0.0276215,-1.00656e-05,1.75773e-09,-1.18886e-13,12491.8,-42.9164], Tmin=(1028.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]CC([CH2])C(=C)[O](2504)',
    structure = SMILES('[CH2]CC([CH2])C(=C)[O]'),
    E0 = (261.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,405.999,406.002,406.003],'cm^-1')),
        HinderedRotor(inertia=(0.006215,'amu*angstrom^2'), symmetry=1, barrier=(0.726941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00621715,'amu*angstrom^2'), symmetry=1, barrier=(0.727194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00621712,'amu*angstrom^2'), symmetry=1, barrier=(0.727284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0582175,'amu*angstrom^2'), symmetry=1, barrier=(6.80968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667782,0.0656988,-4.36726e-05,6.65012e-09,3.58983e-12,31580.5,31.408], Tmin=(100,'K'), Tmax=(953.375,'K')), NASAPolynomial(coeffs=[14.2976,0.0260055,-8.74214e-06,1.46918e-09,-9.8094e-14,28186.7,-37.8661], Tmin=(953.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C(=C)O(4301)',
    structure = SMILES('[CH2]C(O)=C([CH2])[CH]C'),
    E0 = (145.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310478,0.0675362,-2.89799e-05,-2.06124e-08,1.58072e-11,17629,25.8042], Tmin=(100,'K'), Tmax=(940.527,'K')), NASAPolynomial(coeffs=[19.6704,0.0189877,-5.43921e-06,8.97759e-10,-6.32959e-14,12492.8,-74.354], Tmin=(940.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(Allyl_S) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(O)C([CH2])[CH]C(4302)',
    structure = SMILES('[CH]=C(O)C([CH2])[CH]C'),
    E0 = (360.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,573.702],'cm^-1')),
        HinderedRotor(inertia=(0.0101866,'amu*angstrom^2'), symmetry=1, barrier=(2.3604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0569651,'amu*angstrom^2'), symmetry=1, barrier=(13.2185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101398,'amu*angstrom^2'), symmetry=1, barrier=(2.35268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.57358,'amu*angstrom^2'), symmetry=1, barrier=(13.188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.571967,'amu*angstrom^2'), symmetry=1, barrier=(13.207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.385016,0.0711967,-5.35711e-05,1.26026e-08,2.37796e-12,43448.6,32.0949], Tmin=(100,'K'), Tmax=(955.847,'K')), NASAPolynomial(coeffs=[16.4467,0.0229375,-7.5848e-06,1.27584e-09,-8.589e-14,39512.2,-49.198], Tmin=(955.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])CC(2505)',
    structure = SMILES('[CH]=C([O])C([CH2])CC'),
    E0 = (303.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,354.649,354.75],'cm^-1')),
        HinderedRotor(inertia=(0.00134255,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00134224,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110122,'amu*angstrom^2'), symmetry=1, barrier=(9.8247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109894,'amu*angstrom^2'), symmetry=1, barrier=(9.82645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558151,0.0684699,-4.98061e-05,1.16713e-08,2.16423e-12,36617.5,30.4183], Tmin=(100,'K'), Tmax=(947.725,'K')), NASAPolynomial(coeffs=[14.757,0.0254458,-8.46472e-06,1.4103e-09,-9.36142e-14,33167,-41.3257], Tmin=(947.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(C)[CH]C(2744)',
    structure = SMILES('[CH]=C([O])C(C)[CH]C'),
    E0 = (292.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00838956,'amu*angstrom^2'), symmetry=1, barrier=(5.05059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660118,'amu*angstrom^2'), symmetry=1, barrier=(15.1774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.660145,'amu*angstrom^2'), symmetry=1, barrier=(15.178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252088,'amu*angstrom^2'), symmetry=1, barrier=(15.1788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598771,0.0670083,-5.31394e-05,2.21779e-08,-3.74828e-12,35348.2,31.0538], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[14.7529,0.0267181,-1.01315e-05,1.77399e-09,-1.18255e-13,31370.3,-42.0467], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(=C)O(4221)',
    structure = SMILES('[CH2][CH]C([CH2])C(=C)O'),
    E0 = (318.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,886.59],'cm^-1')),
        HinderedRotor(inertia=(0.0755603,'amu*angstrom^2'), symmetry=1, barrier=(1.73728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0755616,'amu*angstrom^2'), symmetry=1, barrier=(1.73731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00311458,'amu*angstrom^2'), symmetry=1, barrier=(1.73726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700329,'amu*angstrom^2'), symmetry=1, barrier=(16.102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288679,'amu*angstrom^2'), symmetry=1, barrier=(16.1019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.496011,0.0684083,-4.73716e-05,7.48712e-09,3.84798e-12,38411.6,33.0798], Tmin=(100,'K'), Tmax=(960.186,'K')), NASAPolynomial(coeffs=[15.9848,0.0235017,-7.86494e-06,1.33537e-09,-9.04256e-14,34532.8,-45.7245], Tmin=(960.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C[C]=O(2467)',
    structure = SMILES('[CH2]C([CH]C)C[C]=O'),
    E0 = (279.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,413.61,2621.97],'cm^-1')),
        HinderedRotor(inertia=(0.0627707,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627708,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627705,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627707,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141566,'amu*angstrom^2'), symmetry=1, barrier=(69.0623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05704,0.0686662,-7.37601e-05,5.27733e-08,-1.65495e-11,33685.5,31.5724], Tmin=(100,'K'), Tmax=(762.453,'K')), NASAPolynomial(coeffs=[6.46712,0.0401417,-1.73637e-05,3.21791e-09,-2.20793e-13,32864.7,6.96619], Tmin=(762.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=O)[CH]C(4303)',
    structure = SMILES('[CH2]C([C]=O)[CH]C'),
    E0 = (312.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,343.549],'cm^-1')),
        HinderedRotor(inertia=(0.00143191,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00142255,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0805419,'amu*angstrom^2'), symmetry=1, barrier=(6.77155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0809389,'amu*angstrom^2'), symmetry=1, barrier=(6.7757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61655,0.0534098,-4.66372e-05,2.28478e-08,-4.73055e-12,37664.9,26.3049], Tmin=(100,'K'), Tmax=(1128.27,'K')), NASAPolynomial(coeffs=[9.03737,0.0271015,-1.16615e-05,2.18172e-09,-1.51458e-13,35990.3,-10.3919], Tmin=(1128.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1C(=O)CC1C(4172)',
    structure = SMILES('[CH2]C1C(=O)CC1C'),
    E0 = (32.2187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36564,0.045488,9.11038e-06,-3.93571e-08,1.68516e-11,3980.72,23.9509], Tmin=(100,'K'), Tmax=(1017.45,'K')), NASAPolynomial(coeffs=[12.1834,0.0315135,-1.23843e-05,2.31013e-09,-1.64029e-13,301.416,-35.689], Tmin=(1017.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.2187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C[CH]C1CCC1=O(3815)',
    structure = SMILES('C[CH]C1CCC1=O'),
    E0 = (30.6516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79069,0.0333136,4.13679e-05,-7.12334e-08,2.7762e-11,3779.57,25.5429], Tmin=(100,'K'), Tmax=(995.412,'K')), NASAPolynomial(coeffs=[11.5285,0.0315103,-1.21636e-05,2.29093e-09,-1.65257e-13,-8.33494,-30.6806], Tmin=(995.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.6516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]C(=CC)C(C)=O(4304)',
    structure = SMILES('[CH2]C(=CC)C(C)=O'),
    E0 = (-49.5112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40859,0.0567725,-3.34456e-05,9.31708e-09,-1.04526e-12,-5862.37,23.3674], Tmin=(100,'K'), Tmax=(1959.15,'K')), NASAPolynomial(coeffs=[14.5088,0.0300261,-1.29678e-05,2.34895e-09,-1.56101e-13,-10995.5,-48.644], Tmin=(1959.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.5112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + radical(C=C(C=O)CJ)"""),
)

species(
    label = '[CH2]C(C=C)C(C)=O(3820)',
    structure = SMILES('[CH2]C(C=C)C(C)=O'),
    E0 = (25.6432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0457,0.0698612,-6.91281e-05,4.25146e-08,-1.15034e-11,3186.32,26.1135], Tmin=(100,'K'), Tmax=(865.523,'K')), NASAPolynomial(coeffs=[7.3609,0.0406752,-1.85464e-05,3.55366e-09,-2.49653e-13,2093.14,-3.44156], Tmin=(865.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.6432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][C]1OC(C)C1[CH2](4305)',
    structure = SMILES('[CH2][C]1OC(C)C1[CH2]'),
    E0 = (394.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.966215,0.0651225,-4.10296e-05,-5.03185e-09,1.3953e-11,47535.6,24.9403], Tmin=(100,'K'), Tmax=(751.294,'K')), NASAPolynomial(coeffs=[11.8592,0.0275596,-6.82866e-06,8.1977e-10,-4.00646e-14,45322.2,-28.3347], Tmin=(751.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(Isobutyl) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1OCC1[CH]C(4306)',
    structure = SMILES('[CH2][C]1OCC1[CH]C'),
    E0 = (396.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29292,0.0623229,-4.18834e-05,1.89214e-10,1.24087e-11,47772.5,25.5972], Tmin=(100,'K'), Tmax=(666.677,'K')), NASAPolynomial(coeffs=[8.20553,0.034929,-1.19305e-05,1.91866e-09,-1.20404e-13,46537.8,-7.2964], Tmin=(666.677,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(Cs_S) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]([CH]C)C(C)=O(4307)',
    structure = SMILES('[CH2]C([CH]C)=C(C)[O]'),
    E0 = (124.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.672031,0.0634948,-3.49297e-05,-6.09796e-10,5.05083e-12,15072.8,25.3182], Tmin=(100,'K'), Tmax=(1016.86,'K')), NASAPolynomial(coeffs=[14.6971,0.0275626,-1.03036e-05,1.85017e-09,-1.28103e-13,11225.9,-47.4699], Tmin=(1016.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)=O(4230)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)=O'),
    E0 = (304.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.437306,'amu*angstrom^2'), symmetry=1, barrier=(10.0545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182105,'amu*angstrom^2'), symmetry=1, barrier=(10.0549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00511079,'amu*angstrom^2'), symmetry=1, barrier=(33.0519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0499496,'amu*angstrom^2'), symmetry=1, barrier=(33.0523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43754,'amu*angstrom^2'), symmetry=1, barrier=(33.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02005,0.0696661,-6.71767e-05,3.87794e-08,-9.69296e-12,36683.9,30.6792], Tmin=(100,'K'), Tmax=(937.74,'K')), NASAPolynomial(coeffs=[8.33828,0.0384505,-1.72459e-05,3.28308e-09,-2.29969e-13,35311.4,-4.15681], Tmin=(937.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(CJC(C)C=O) + radical(RCCJ)"""),
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
    E0 = (250.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (393.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (408.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (603.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (608.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (720.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1084.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (716.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (705.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (258.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (259.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (259.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (273.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (273.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (314.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (275.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (275.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (650.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (482.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (387.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (386.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (408.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (406.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (318.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (397.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (378.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (377.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (303.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (435.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (611.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (494.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (667.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (709.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (387.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (392.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (413.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (407.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (552.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (347.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (337.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (391.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (524.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (750.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (259.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (259.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (314.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (275.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (395.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (396.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (367.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (358.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['ketene(1375)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=C([O])C[CH][CH]C(2560)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C(=C)[O](2379)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', 'C=C([O])[CH][CH]C(2452)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', '[CH]C([CH2])C(=C)[O](2593)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH2]C([C]=C)[CH]C(4287)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2]C([C]C)C(=C)[O](4288)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([CH]C)C(=C)[O](4289)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=C([O])C1CC1C(4290)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C1C(=C)OC1C(4188)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=C1OCC1[CH]C(4273)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=C([O])C(=C)CC(2501)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=C([O])C(C)=CC(2740)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C(=CC)C(=C)O(4291)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=CC(C)C(=C)[O](2741)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C(C=C)C(=C)O(3819)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]([CH]C)C([CH2])[O](4292)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C([CH]C)[C]1CO1(4293)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C1[C]([O])CC1C(4294)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.16207e+08,'s^-1'), n=0.911389, Ea=(137.014,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 135.1 to 137.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C[CH]C1CC[C]1[O](4295)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(135.515,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 133.4 to 135.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C1C(C)C1([CH2])[O](4212)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(157.645,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C1([O])CC1[CH]C(4296)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(155.929,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C(=CC)C(=C)[O](4297)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]C(C=C)C(=C)[O](4298)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(T)(20)', 'C=C([O])C=CC(2451)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.898708,'m^3/(mol*s)'), n=2.12172, Ea=(14.1789,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-CsH;YJ] for rate rule [Cds-CdH_Cds-CsH;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CHCH3(T)(21)', 'C=CC(=C)[O](2319)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.0532,'m^3/(mol*s)'), n=2.065, Ea=(15.9938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;YJ] for rate rule [Cds-CdH_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=O(1376)', 'm1_allyl(186)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['ketene(1375)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]=O(1376)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][C]([CH]C)C(=C)[O](4299)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C(=C)[O](4216)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH]=C([O])C([CH2])[CH]C(4300)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2][C](CC)C(=C)[O](2502)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C=C([O])[C](C)[CH]C(2742)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC([CH2])C(=C)[O](2504)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2][C]([CH]C)C(=C)O(4301)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C(O)C([CH2])[CH]C(4302)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C([O])C([CH2])CC(2505)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C([O])C(C)[CH]C(2744)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])C(=C)O(4221)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CH2(T)(20)', '[CH2]C([C]=O)[CH]C(4303)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C1C(=O)CC1C(4172)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['C[CH]C1CCC1=O(3815)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C(=CC)C(C)=O(4304)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C(C=C)C(C)=O(3820)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2][C]1OC(C)C1[CH2](4305)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(144.316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2][C]1OCC1[CH]C(4306)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(146.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2][C]([CH]C)C(C)=O(4307)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][CH]C([CH2])C(C)=O(4230)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1193',
    isomers = [
        '[CH2]C([CH]C)C(=C)[O](2503)',
    ],
    reactants = [
        ('ketene(1375)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1193',
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


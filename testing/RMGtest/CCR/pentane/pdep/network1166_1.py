species(
    label = '[CH2]C([O])C([CH2])[CH]C(1769)',
    structure = SMILES('[CH2]C([O])C([CH2])[CH]C'),
    E0 = (474.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200.373,800.746,1602.98],'cm^-1')),
        HinderedRotor(inertia=(0.155436,'amu*angstrom^2'), symmetry=1, barrier=(3.58712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155436,'amu*angstrom^2'), symmetry=1, barrier=(3.58712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155436,'amu*angstrom^2'), symmetry=1, barrier=(3.58712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155436,'amu*angstrom^2'), symmetry=1, barrier=(3.58712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155436,'amu*angstrom^2'), symmetry=1, barrier=(3.58712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679807,0.0728465,-6.33018e-05,3.12085e-08,-6.472e-12,57236,32.8704], Tmin=(100,'K'), Tmax=(1133.4,'K')), NASAPolynomial(coeffs=[11.054,0.0362334,-1.48458e-05,2.70638e-09,-1.85091e-13,54884.4,-18.4783], Tmin=(1133.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3817.96,'J/mol'), sigma=(6.83058,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.36 K, Pc=27.18 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.824856,0.0704451,-5.62754e-05,2.453e-08,-4.52532e-12,57249.7,32.3517], Tmin=(100,'K'), Tmax=(1246.76,'K')), NASAPolynomial(coeffs=[11.2555,0.0369797,-1.6012e-05,3.00011e-09,-2.08094e-13,54648.8,-20.2706], Tmin=(1246.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([O])C[CH][CH]C(1815)',
    structure = SMILES('[CH2]C([O])C[CH][CH]C'),
    E0 = (466.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,359.638,441.319,1889.65,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83996,0.0764746,-9.64501e-05,8.3381e-08,-3.00819e-11,56216.4,33.5352], Tmin=(100,'K'), Tmax=(821.465,'K')), NASAPolynomial(coeffs=[4.10855,0.0480153,-2.15791e-05,4.03073e-09,-2.75913e-13,56102.6,20.9849], Tmin=(821.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.89,'J/mol'), sigma=(6.79962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.22 K, Pc=27.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.685187,0.0807927,-0.000107586,9.46849e-08,-3.41163e-11,54425.3,33.1627], Tmin=(100,'K'), Tmax=(830.633,'K')), NASAPolynomial(coeffs=[4.28037,0.0484164,-2.19173e-05,4.09419e-09,-2.79718e-13,54347.7,19.6132], Tmin=(830.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH)"""),
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
    label = '[CH2][CH]C([CH2])[O](1582)',
    structure = SMILES('[CH2][CH]C([CH2])[O]'),
    E0 = (530.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1127.94],'cm^-1')),
        HinderedRotor(inertia=(0.109106,'amu*angstrom^2'), symmetry=1, barrier=(2.50857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0027686,'amu*angstrom^2'), symmetry=1, barrier=(2.49694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0182808,'amu*angstrom^2'), symmetry=1, barrier=(16.5026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78639,0.0431654,-3.16516e-05,9.87645e-09,-5.24257e-13,63867.6,24.633], Tmin=(100,'K'), Tmax=(1082.15,'K')), NASAPolynomial(coeffs=[11.0036,0.0170489,-6.47571e-06,1.15861e-09,-7.93421e-14,61407,-22.7143], Tmin=(1082.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CJCO) + radical(RCCJ)"""),
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
    label = '[CH2]C([O])[CH][CH]C(1681)',
    structure = SMILES('[CH2]C([O])[CH][CH]C'),
    E0 = (495.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,368.964,1232.71,1232.73],'cm^-1')),
        HinderedRotor(inertia=(0.0632841,'amu*angstrom^2'), symmetry=1, barrier=(6.11349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0632871,'amu*angstrom^2'), symmetry=1, barrier=(6.11351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00123841,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00566925,'amu*angstrom^2'), symmetry=1, barrier=(6.11344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57704,0.0541882,-4.32836e-05,1.99115e-08,-3.9589e-12,59710.6,28.5786], Tmin=(100,'K'), Tmax=(1154.82,'K')), NASAPolynomial(coeffs=[8.32727,0.0308074,-1.29145e-05,2.37992e-09,-1.63627e-13,58151.5,-4.95905], Tmin=(1154.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(CJCO)"""),
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
    label = '[CH2][CH]C([CH2])[CH]C(3871)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]C'),
    E0 = (605.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,246.588,3146],'cm^-1')),
        HinderedRotor(inertia=(0.00277254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.6228e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0110378,'amu*angstrom^2'), symmetry=1, barrier=(28.5633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00283321,'amu*angstrom^2'), symmetry=1, barrier=(7.3317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07418,'amu*angstrom^2'), symmetry=1, barrier=(46.3499,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90503,0.0506182,-2.74878e-05,7.53461e-09,-8.92238e-13,72859.4,29.9885], Tmin=(100,'K'), Tmax=(1669.08,'K')), NASAPolynomial(coeffs=[7.3347,0.0376059,-1.57937e-05,2.86374e-09,-1.92623e-13,71046.9,1.01198], Tmin=(1669.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][O](4436)',
    structure = SMILES('[CH2]C([CH]C)[CH][O]'),
    E0 = (472.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1413.87,1414.17,1414.31],'cm^-1')),
        HinderedRotor(inertia=(0.179074,'amu*angstrom^2'), symmetry=1, barrier=(4.11726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178866,'amu*angstrom^2'), symmetry=1, barrier=(4.11247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17914,'amu*angstrom^2'), symmetry=1, barrier=(4.11877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179155,'amu*angstrom^2'), symmetry=1, barrier=(4.11913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36092,0.0656478,-9.23383e-05,8.44737e-08,-3.08504e-11,56858.8,28.4924], Tmin=(100,'K'), Tmax=(848.262,'K')), NASAPolynomial(coeffs=[3.29806,0.0401068,-1.81616e-05,3.3756e-09,-2.29214e-13,57120.4,22.9448], Tmin=(848.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(472.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(CCsJOH) + radical(Isobutyl)"""),
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
    label = '[CH]C([CH2])C([CH2])[O](1869)',
    structure = SMILES('[CH]C([CH2])C([CH2])[O]'),
    E0 = (748.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8706,0.0648702,-6.65585e-05,3.64488e-08,-7.92718e-12,90199.8,28.3856], Tmin=(100,'K'), Tmax=(1122.48,'K')), NASAPolynomial(coeffs=[13.3301,0.0204705,-7.22632e-06,1.21023e-09,-7.88509e-14,87402.6,-33.164], Tmin=(1122.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(748.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C([O])C([CH2])[C]C(4437)',
    structure = SMILES('[CH2]C([O])C([CH2])[C]C'),
    E0 = (728.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.261411,0.079211,-7.88582e-05,4.2528e-08,-9.23509e-12,87763.4,31.7194], Tmin=(100,'K'), Tmax=(1114.75,'K')), NASAPolynomial(coeffs=[14.3228,0.0287546,-1.09638e-05,1.92387e-09,-1.28902e-13,84628.4,-37.6462], Tmin=(1114.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(728.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH]C)C([CH2])[O](4438)',
    structure = SMILES('[CH]C([CH]C)C([CH2])[O]'),
    E0 = (718.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721053,0.0728902,-6.56319e-05,3.2254e-08,-6.60016e-12,86475.9,31.9456], Tmin=(100,'K'), Tmax=(1149.87,'K')), NASAPolynomial(coeffs=[11.9189,0.0339366,-1.48169e-05,2.79244e-09,-1.94711e-13,83900.7,-23.641], Tmin=(1149.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(718.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])C([CH2])[CH]C(4439)',
    structure = SMILES('[CH]C([O])C([CH2])[CH]C'),
    E0 = (711.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.78045,0.0706981,-6.29954e-05,3.16117e-08,-6.63804e-12,85691.9,32.8799], Tmin=(100,'K'), Tmax=(1123.58,'K')), NASAPolynomial(coeffs=[11.0793,0.034034,-1.40483e-05,2.56942e-09,-1.76072e-13,83377.6,-18.006], Tmin=(1123.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(711.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])C1CC1C(4440)',
    structure = SMILES('[CH2]C([O])C1CC1C'),
    E0 = (222.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894893,0.0512386,1.58225e-05,-5.96652e-08,2.71867e-11,26843,26.3381], Tmin=(100,'K'), Tmax=(966.554,'K')), NASAPolynomial(coeffs=[17.0565,0.0257307,-8.80271e-06,1.6083e-09,-1.17041e-13,21786.1,-61.0807], Tmin=(966.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH]C)C1CO1(4441)',
    structure = SMILES('[CH2]C([CH]C)C1CO1'),
    E0 = (219.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19606,0.0679129,-5.14126e-05,2.25556e-08,-3.94882e-12,26509.2,31.2094], Tmin=(100,'K'), Tmax=(1587.07,'K')), NASAPolynomial(coeffs=[12.3568,0.0293413,-7.46963e-06,9.51705e-10,-5.02673e-14,23646.9,-29.9328], Tmin=(1587.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1OC(C)C1[CH2](4337)',
    structure = SMILES('[CH2]C1OC(C)C1[CH2]'),
    E0 = (213.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91879,0.0564676,5.95056e-07,-5.11724e-08,2.89233e-11,25812.6,24.452], Tmin=(100,'K'), Tmax=(856.862,'K')), NASAPolynomial(coeffs=[15.2739,0.0244146,-4.49291e-06,4.01238e-10,-1.62062e-14,22069.2,-50.0745], Tmin=(856.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(CJC(C)OC)"""),
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
    label = '[CH2]C1OCC1[CH]C(3810)',
    structure = SMILES('[CH2]C1OCC1[CH]C'),
    E0 = (215.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07239,0.0561098,-1.11714e-05,-2.7065e-08,1.63162e-11,26056.8,25.7071], Tmin=(100,'K'), Tmax=(889.226,'K')), NASAPolynomial(coeffs=[11.7997,0.0314563,-9.39591e-06,1.45137e-09,-9.23885e-14,23215.9,-30.0333], Tmin=(889.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Cs_S) + radical(CJC(C)OC)"""),
)

species(
    label = 'C[CH]C1CCC1[O](4349)',
    structure = SMILES('C[CH]C1CCC1[O]'),
    E0 = (209.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47866,0.0387094,3.91939e-05,-7.27288e-08,2.8699e-11,25327,26.4613], Tmin=(100,'K'), Tmax=(1001.41,'K')), NASAPolynomial(coeffs=[13.0875,0.0329331,-1.29579e-05,2.46861e-09,-1.79156e-13,20966.6,-39.7235], Tmin=(1001.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C([O])C(=C)CC(1765)',
    structure = SMILES('[CH2]C([O])C(=C)CC'),
    E0 = (190.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,756.217],'cm^-1')),
        HinderedRotor(inertia=(0.768498,'amu*angstrom^2'), symmetry=1, barrier=(17.6693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130631,'amu*angstrom^2'), symmetry=1, barrier=(17.6691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0435255,'amu*angstrom^2'), symmetry=1, barrier=(17.6686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0015562,'amu*angstrom^2'), symmetry=1, barrier=(17.6691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.3965,0.0687613,-4.15728e-05,4.78089e-09,3.04214e-12,23101.2,30.3923], Tmin=(100,'K'), Tmax=(1063.33,'K')), NASAPolynomial(coeffs=[15.8668,0.0289465,-1.13369e-05,2.08092e-09,-1.45217e-13,18772.1,-50.0791], Tmin=(1063.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])C(C)=CC(2026)',
    structure = SMILES('[CH2]C([O])C(C)=CC'),
    E0 = (177.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,329.719,329.783],'cm^-1')),
        HinderedRotor(inertia=(0.00155045,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141875,'amu*angstrom^2'), symmetry=1, barrier=(10.9487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141901,'amu*angstrom^2'), symmetry=1, barrier=(10.9491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141728,'amu*angstrom^2'), symmetry=1, barrier=(10.9434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506095,0.0711669,-5.49283e-05,2.21037e-08,-3.63545e-12,21483.5,28.8939], Tmin=(100,'K'), Tmax=(1423.35,'K')), NASAPolynomial(coeffs=[14.7887,0.0310294,-1.26299e-05,2.29236e-09,-1.55798e-13,17417.6,-45.0537], Tmin=(1423.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH]C)C(=C)O(2509)',
    structure = SMILES('[CH2]C([CH]C)C(=C)O'),
    E0 = (113.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,720.45,720.459],'cm^-1')),
        HinderedRotor(inertia=(0.0325029,'amu*angstrom^2'), symmetry=1, barrier=(11.9738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138494,'amu*angstrom^2'), symmetry=1, barrier=(3.18426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325085,'amu*angstrom^2'), symmetry=1, barrier=(11.9737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325026,'amu*angstrom^2'), symmetry=1, barrier=(11.9738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.520769,'amu*angstrom^2'), symmetry=1, barrier=(11.9735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456001,0.067707,-3.70287e-05,-4.63309e-09,8.23313e-12,13729,31.4369], Tmin=(100,'K'), Tmax=(958.982,'K')), NASAPolynomial(coeffs=[16.2161,0.0257542,-8.60957e-06,1.47243e-09,-1.00578e-13,9612.63,-49.6384], Tmin=(958.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)C(C)=O(2519)',
    structure = SMILES('[CH2]C([CH]C)C(C)=O'),
    E0 = (98.8978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,317.875,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00166833,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110587,'amu*angstrom^2'), symmetry=1, barrier=(7.92951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110587,'amu*angstrom^2'), symmetry=1, barrier=(7.9295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339854,'amu*angstrom^2'), symmetry=1, barrier=(24.3688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114749,'amu*angstrom^2'), symmetry=1, barrier=(7.92951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01231,0.0687032,-5.64122e-05,2.66687e-08,-5.49408e-12,11999.7,28.9119], Tmin=(100,'K'), Tmax=(1109.56,'K')), NASAPolynomial(coeffs=[8.98899,0.0399463,-1.7535e-05,3.30908e-09,-2.30689e-13,10229.6,-10.4], Tmin=(1109.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.8978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(=CC)C([CH2])O(4442)',
    structure = SMILES('[CH2]C(=CC)C([CH2])O'),
    E0 = (98.6792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0724147,0.0757481,-6.26083e-05,2.66347e-08,-4.52234e-12,12018.6,30.9305], Tmin=(100,'K'), Tmax=(1413.95,'K')), NASAPolynomial(coeffs=[17.6959,0.0258918,-9.71779e-06,1.69713e-09,-1.13123e-13,7034.89,-60.1974], Tmin=(1413.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.6792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=CC)C(C)[O](4443)',
    structure = SMILES('[CH2]C(=CC)C(C)[O]'),
    E0 = (117.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558841,0.0646756,-3.09708e-05,-4.41007e-09,5.77488e-12,14259.4,28.3506], Tmin=(100,'K'), Tmax=(1059.75,'K')), NASAPolynomial(coeffs=[15.0773,0.0306369,-1.21771e-05,2.2532e-09,-1.5797e-13,10016.4,-48.0356], Tmin=(1059.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=O)C([CH2])CC(1766)',
    structure = SMILES('[CH2]C(CC)C(=C)[O]'),
    E0 = (56.2744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,393.74,393.948,394.399],'cm^-1')),
        HinderedRotor(inertia=(0.00108189,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0010854,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0827327,'amu*angstrom^2'), symmetry=1, barrier=(9.13182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828341,'amu*angstrom^2'), symmetry=1, barrier=(9.13032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.58,'J/mol'), sigma=(6.61681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.14 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626666,0.0650109,-3.33787e-05,-5.40409e-09,7.94584e-12,6897.98,29.7691], Tmin=(100,'K'), Tmax=(953.666,'K')), NASAPolynomial(coeffs=[14.5328,0.0282512,-9.48282e-06,1.60529e-09,-1.08168e-13,3264.88,-41.8022], Tmin=(953.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.2744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=O)C(C)[CH]C(2027)',
    structure = SMILES('C=C([O])C(C)[CH]C'),
    E0 = (45.7342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,180,779.905,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334751,'amu*angstrom^2'), symmetry=1, barrier=(14.4511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00848413,'amu*angstrom^2'), symmetry=1, barrier=(3.66081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628529,'amu*angstrom^2'), symmetry=1, barrier=(14.4511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62853,'amu*angstrom^2'), symmetry=1, barrier=(14.4511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.58,'J/mol'), sigma=(6.61681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.14 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722002,0.0629666,-3.49958e-05,3.37878e-09,2.54443e-12,5626.3,30.2046], Tmin=(100,'K'), Tmax=(1084.24,'K')), NASAPolynomial(coeffs=[13.4125,0.0312882,-1.2115e-05,2.18829e-09,-1.50469e-13,1984.5,-36.1499], Tmin=(1084.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C([O])C(C)C=C(2028)',
    structure = SMILES('[CH2]C([O])C(C)C=C'),
    E0 = (197.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,477.782,1197.87],'cm^-1')),
        HinderedRotor(inertia=(0.00173224,'amu*angstrom^2'), symmetry=1, barrier=(19.668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855408,'amu*angstrom^2'), symmetry=1, barrier=(19.6675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85545,'amu*angstrom^2'), symmetry=1, barrier=(19.6685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855424,'amu*angstrom^2'), symmetry=1, barrier=(19.6679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449935,0.0681861,-4.18444e-05,6.25139e-09,2.34802e-12,23881,29.8658], Tmin=(100,'K'), Tmax=(1066.23,'K')), NASAPolynomial(coeffs=[15.1528,0.0298395,-1.15484e-05,2.09653e-09,-1.45092e-13,19790.1,-46.4912], Tmin=(1066.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C([CH2])C=C(3812)',
    structure = SMILES('[CH2]C(O)C([CH2])C=C'),
    E0 = (172.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.167345,0.0754739,-6.57891e-05,3.08102e-08,-5.78923e-12,20849.9,32.8542], Tmin=(100,'K'), Tmax=(1285.72,'K')), NASAPolynomial(coeffs=[15.712,0.0271128,-9.36777e-06,1.55464e-09,-1.00661e-13,16852.7,-46.0463], Tmin=(1285.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(C=C)C(C)[O](3814)',
    structure = SMILES('[CH2]C(C=C)C(C)[O]'),
    E0 = (190.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.657187,0.0642814,-3.33716e-05,-1.72822e-09,5.33503e-12,23090.6,30.2677], Tmin=(100,'K'), Tmax=(1007.6,'K')), NASAPolynomial(coeffs=[13.7696,0.0308376,-1.12895e-05,1.99216e-09,-1.36184e-13,19503.5,-37.7795], Tmin=(1007.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=CC)C([CH2])[O](4444)',
    structure = SMILES('[CH2]C(=CC)C([CH2])[O]'),
    E0 = (329.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,277.926],'cm^-1')),
        HinderedRotor(inertia=(0.00401278,'amu*angstrom^2'), symmetry=1, barrier=(2.40336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768385,'amu*angstrom^2'), symmetry=1, barrier=(17.6667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76799,'amu*angstrom^2'), symmetry=1, barrier=(17.6576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0294772,'amu*angstrom^2'), symmetry=1, barrier=(17.6664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.365304,0.0703335,-5.2269e-05,1.66162e-08,-1.11365e-12,39713.3,29.6801], Tmin=(100,'K'), Tmax=(1114.01,'K')), NASAPolynomial(coeffs=[16.1801,0.0263861,-1.03801e-05,1.89265e-09,-1.30891e-13,35393.2,-51.8999], Tmin=(1114.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.58,'J/mol'), sigma=(6.61681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.14 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586634,0.0664197,-5.34853e-05,2.33575e-08,-4.13471e-12,30296.4,32.4745], Tmin=(100,'K'), Tmax=(1352.09,'K')), NASAPolynomial(coeffs=[14.0751,0.0265156,-9.2158e-06,1.52979e-09,-9.87714e-14,26648.9,-36.6682], Tmin=(1352.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([O])C([CH2])C=C(4445)',
    structure = SMILES('[CH2]C([O])C([CH2])C=C'),
    E0 = (402.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,604.318,3113.3],'cm^-1')),
        HinderedRotor(inertia=(0.0720765,'amu*angstrom^2'), symmetry=1, barrier=(18.6837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.52797,'amu*angstrom^2'), symmetry=1, barrier=(81.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721023,'amu*angstrom^2'), symmetry=1, barrier=(18.6838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139872,'amu*angstrom^2'), symmetry=1, barrier=(3.21594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491081,0.0696401,-5.37465e-05,1.82912e-08,-1.21122e-12,48543.4,31.4972], Tmin=(100,'K'), Tmax=(1021.56,'K')), NASAPolynomial(coeffs=[14.6069,0.0270162,-9.73121e-06,1.68648e-09,-1.13564e-13,44999.4,-40.1346], Tmin=(1021.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])C=CC(1678)',
    structure = SMILES('[CH2]C([O])C=CC'),
    E0 = (216.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,235.724,4000],'cm^-1')),
        HinderedRotor(inertia=(0.487405,'amu*angstrom^2'), symmetry=1, barrier=(18.8908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488885,'amu*angstrom^2'), symmetry=1, barrier=(18.8911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.480351,'amu*angstrom^2'), symmetry=1, barrier=(18.8731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23576,0.0525332,-2.93684e-05,1.95313e-09,2.57797e-12,26157,25.5176], Tmin=(100,'K'), Tmax=(1087.33,'K')), NASAPolynomial(coeffs=[12.793,0.0241196,-9.62576e-06,1.77652e-09,-1.23939e-13,22810,-35.0406], Tmin=(1087.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])C=C(1523)',
    structure = SMILES('[CH2]C([O])C=C'),
    E0 = (252.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,529.328,529.355],'cm^-1')),
        HinderedRotor(inertia=(0.0692523,'amu*angstrom^2'), symmetry=1, barrier=(13.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0692738,'amu*angstrom^2'), symmetry=1, barrier=(13.7796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95358,0.036167,-6.30157e-06,-1.92918e-08,1.03178e-11,30464.9,21.1322], Tmin=(100,'K'), Tmax=(988.725,'K')), NASAPolynomial(coeffs=[12.0826,0.0154802,-5.70142e-06,1.06015e-09,-7.65757e-14,27470.1,-32.6352], Tmin=(988.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH]C)C=O(3827)',
    structure = SMILES('[CH2]C([CH]C)C=O'),
    E0 = (153.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,298.571],'cm^-1')),
        HinderedRotor(inertia=(0.093393,'amu*angstrom^2'), symmetry=1, barrier=(5.91021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00189111,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00190494,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947027,'amu*angstrom^2'), symmetry=1, barrier=(5.90918,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62239,0.0540134,-4.22431e-05,1.86389e-08,-3.56669e-12,18578.2,25.1257], Tmin=(100,'K'), Tmax=(1186.46,'K')), NASAPolynomial(coeffs=[8.33857,0.0313704,-1.36159e-05,2.55323e-09,-1.77214e-13,16984.5,-8.4242], Tmin=(1186.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C([CH]C)C=C(3803)',
    structure = SMILES('[CH2]C([CH]C)C=C'),
    E0 = (327.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1026.84,3543.87],'cm^-1')),
        HinderedRotor(inertia=(0.235199,'amu*angstrom^2'), symmetry=1, barrier=(5.40768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235218,'amu*angstrom^2'), symmetry=1, barrier=(5.40812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271982,'amu*angstrom^2'), symmetry=1, barrier=(20.3613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.76263,'amu*angstrom^2'), symmetry=1, barrier=(86.5103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39619,0.0496388,-1.64699e-05,-7.33995e-09,4.82136e-12,39497.8,28.008], Tmin=(100,'K'), Tmax=(1082.86,'K')), NASAPolynomial(coeffs=[9.7195,0.0328979,-1.26799e-05,2.27041e-09,-1.54839e-13,36874.1,-16.601], Tmin=(1082.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH2][C]([O])C([CH2])[CH]C(4446)',
    structure = SMILES('[CH2][C]([O])C([CH2])[CH]C'),
    E0 = (651.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,245.127,884.613,1961.02],'cm^-1')),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446731,0.0824915,-0.000103553,7.7858e-08,-2.42038e-11,78483.6,33.1335], Tmin=(100,'K'), Tmax=(804.105,'K')), NASAPolynomial(coeffs=[9.30196,0.0369728,-1.59014e-05,2.91671e-09,-1.98086e-13,77107,-7.36197], Tmin=(804.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(Isobutyl) + radical(CJCO)"""),
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
    label = '[CH2][C](CC)C([CH2])[O](1768)',
    structure = SMILES('[CH2][C](CC)C([CH2])[O]'),
    E0 = (432.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.348222,0.0722288,-5.71168e-05,2.43114e-08,-4.22027e-12,52206.9,31.7016], Tmin=(100,'K'), Tmax=(1367.31,'K')), NASAPolynomial(coeffs=[14.6552,0.0303736,-1.11991e-05,1.92272e-09,-1.26631e-13,48294.5,-41.7971], Tmin=(1367.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(CJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([O])[C](C)[CH]C(2030)',
    structure = SMILES('[CH2]C([O])[C](C)[CH]C'),
    E0 = (422.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772642,0.0663494,-4.54507e-05,1.59358e-08,-2.30332e-12,50920.7,30.9523], Tmin=(100,'K'), Tmax=(1575.86,'K')), NASAPolynomial(coeffs=[13.952,0.0328965,-1.36082e-05,2.46491e-09,-1.66251e-13,46766.9,-38.6244], Tmin=(1575.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO)"""),
)

species(
    label = '[CH2][C](O)C([CH2])[CH]C(2515)',
    structure = SMILES('[CH2][C](O)C([CH2])[CH]C'),
    E0 = (421.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,202.936,1767.63],'cm^-1')),
        HinderedRotor(inertia=(0.118936,'amu*angstrom^2'), symmetry=1, barrier=(3.34072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118936,'amu*angstrom^2'), symmetry=1, barrier=(3.34072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118936,'amu*angstrom^2'), symmetry=1, barrier=(3.34072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118936,'amu*angstrom^2'), symmetry=1, barrier=(3.34072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118936,'amu*angstrom^2'), symmetry=1, barrier=(3.34072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118936,'amu*angstrom^2'), symmetry=1, barrier=(3.34072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.265999,0.0867131,-0.00011033,8.41206e-08,-2.63714e-11,50783.9,33.9729], Tmin=(100,'K'), Tmax=(825.989,'K')), NASAPolynomial(coeffs=[9.44488,0.0386276,-1.64051e-05,2.98434e-09,-2.01396e-13,49391.5,-7.80447], Tmin=(825.989,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(C2CsJOH) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH]C)[C](C)[O](2496)',
    structure = SMILES('[CH2]C([CH]C)[C](C)[O]'),
    E0 = (439.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,360,370,350,3000,3100,440,815,1455,1000,206.199,818.639,1674.55],'cm^-1')),
        HinderedRotor(inertia=(0.145708,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145708,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145708,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145708,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145708,'amu*angstrom^2'), symmetry=1, barrier=(3.56099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74109,0.0755035,-7.6901e-05,4.8798e-08,-1.33615e-11,53025.5,31.4514], Tmin=(100,'K'), Tmax=(865.916,'K')), NASAPolynomial(coeffs=[8.12465,0.0413955,-1.78163e-05,3.30833e-09,-2.27943e-13,51746.8,-3.10694], Tmin=(865.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])[O](1771)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])[O]'),
    E0 = (485.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356094,0.0768309,-6.95752e-05,3.4756e-08,-7.10384e-12,58537.8,33.2603], Tmin=(100,'K'), Tmax=(1170.32,'K')), NASAPolynomial(coeffs=[13.4058,0.0322282,-1.24072e-05,2.19022e-09,-1.47157e-13,55483.3,-31.7493], Tmin=(1170.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])C([CH2])CC(1770)',
    structure = SMILES('[CH2][C]([O])C([CH2])CC'),
    E0 = (456.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.31349,0.0829813,-8.93227e-05,5.54369e-08,-1.41824e-11,55092.9,31.0585], Tmin=(100,'K'), Tmax=(942.326,'K')), NASAPolynomial(coeffs=[11.3551,0.0361124,-1.47181e-05,2.65745e-09,-1.80196e-13,53011.9,-21.5554], Tmin=(942.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])C(C)[CH]C(2031)',
    structure = SMILES('[CH2][C]([O])C(C)[CH]C'),
    E0 = (446.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,360,370,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695699,0.0774796,-7.85736e-05,4.79491e-08,-1.2603e-11,53808.9,30.4707], Tmin=(100,'K'), Tmax=(897.79,'K')), NASAPolynomial(coeffs=[8.72156,0.0417229,-1.8835e-05,3.5914e-09,-2.51623e-13,52367.8,-7.3845], Tmin=(897.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([CH]C)C([CH2])O(4447)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])O'),
    E0 = (397.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767967,0.0704535,-5.88535e-05,2.78934e-08,-5.58428e-12,47877.6,32.9397], Tmin=(100,'K'), Tmax=(1167.59,'K')), NASAPolynomial(coeffs=[10.8087,0.0360549,-1.4661e-05,2.66016e-09,-1.81345e-13,45533,-17.0565], Tmin=(1167.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH]C)C(C)[O](4448)',
    structure = SMILES('[CH2][C]([CH]C)C(C)[O]'),
    E0 = (415.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.839031,0.0640612,-4.24494e-05,1.48279e-08,-2.1587e-12,50136.5,31.8629], Tmin=(100,'K'), Tmax=(1559.62,'K')), NASAPolynomial(coeffs=[12.7005,0.0336401,-1.31915e-05,2.32157e-09,-1.54014e-13,46436.6,-30.6338], Tmin=(1559.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])O(4412)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])O'),
    E0 = (449.778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,201.369,2438.34],'cm^-1')),
        HinderedRotor(inertia=(0.0881888,'amu*angstrom^2'), symmetry=1, barrier=(2.44494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881888,'amu*angstrom^2'), symmetry=1, barrier=(2.44494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881888,'amu*angstrom^2'), symmetry=1, barrier=(2.44494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881888,'amu*angstrom^2'), symmetry=1, barrier=(2.44494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881888,'amu*angstrom^2'), symmetry=1, barrier=(2.44494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0881888,'amu*angstrom^2'), symmetry=1, barrier=(2.44494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.581361,0.0771737,-7.79394e-05,4.5924e-08,-1.13138e-11,54217.2,35.2081], Tmin=(100,'K'), Tmax=(969.761,'K')), NASAPolynomial(coeffs=[10.4204,0.0365908,-1.5168e-05,2.77221e-09,-1.89609e-13,52308.9,-11.9578], Tmin=(969.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)[O](4413)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)[O]'),
    E0 = (468.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884591,0.0682369,-5.35431e-05,2.37339e-08,-4.49341e-12,56465.8,33.2858], Tmin=(100,'K'), Tmax=(1216.44,'K')), NASAPolynomial(coeffs=[10.304,0.0372631,-1.53491e-05,2.80179e-09,-1.91489e-13,54174.1,-14.0031], Tmin=(1216.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (474.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (618.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (632.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (632.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (929.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (933.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1124.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (909.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (940.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (940.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (929.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (923.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (482.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (480.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (483.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (483.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (483.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (483.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (497.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (497.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (497.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (497.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (538.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (538.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (538.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (538.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (499.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (499.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (499.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (552.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (484.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (616.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (599.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (596.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (592.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (570.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (510.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (514.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (818.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (839.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (863.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (891.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (611.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (616.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (588.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (616.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (633.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (618.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (592.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (550.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (592.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (534.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (529.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['vinoxy(1351)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C(C)C([CH2])[O](2032)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([O])C[CH][CH]C(1815)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([CH]C)C[CH][O](1700)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C([CH2])[O](1582)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2]C([O])[CH][CH]C(1681)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH2][CH]C([CH2])[CH]C(3871)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH2]C([CH]C)[CH][O](4436)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH]C([CH2])C([CH2])[O](1869)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2]C([O])C([CH2])[C]C(4437)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]C([CH]C)C([CH2])[O](4438)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C([O])C([CH2])[CH]C(4439)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([O])C1CC1C(4440)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([CH]C)C1CO1(4441)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C1OC(C)C1[CH2](4337)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C1C(C)CC1[O](4368)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C1OCC1[CH]C(3810)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['C[CH]C1CCC1[O](4349)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([O])C(=C)CC(1765)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([O])C(C)=CC(2026)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([CH]C)C(=C)O(2509)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([CH]C)C(C)=O(2519)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C(=CC)C([CH2])O(4442)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C(=CC)C(C)[O](4443)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C(=O)C([CH2])CC(1766)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C(=O)C(C)[CH]C(2027)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([O])C(C)C=C(2028)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C(O)C([CH2])C=C(3812)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C(C=C)C(C)[O](3814)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2]C(=CC)C([CH2])[O](4444)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]C([CH]C)C(=C)[O](2503)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]C([O])C([CH2])C=C(4445)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(T)(20)', '[CH2]C([O])C=CC(1678)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(41.7,'m^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;Y_1centerbirad] for rate rule [Cds-CsH_Cds-CsH;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['CHCH3(T)(21)', '[CH2]C([O])C=C(1523)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(T)(20)', '[CH2]C([CH]C)C=O(3827)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(4)', '[CH2]C([CH]C)C=C(3803)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH][O](1556)', 'm1_allyl(186)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['vinoxy(1351)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', '[CH2][C]([CH]C)C([CH2])[O](4292)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(3)', '[CH2][C]([O])C([CH2])[CH]C(4446)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][C](CC)C([CH2])[O](1768)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([O])[C](C)[CH]C(2030)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][C](O)C([CH2])[CH]C(2515)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2]C([CH]C)[C](C)[O](2496)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]CC([CH2])C([CH2])[O](1771)'],
    products = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][C]([O])C([CH2])CC(1770)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][C]([O])C(C)[CH]C(2031)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][C]([CH]C)C([CH2])O(4447)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][C]([CH]C)C(C)[O](4448)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][CH]C([CH2])C([CH2])O(4412)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C([O])C([CH2])[CH]C(1769)'],
    products = ['[CH2][CH]C([CH2])C(C)[O](4413)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1166',
    isomers = [
        '[CH2]C([O])C([CH2])[CH]C(1769)',
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
    label = 'PDepNetwork #1166',
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


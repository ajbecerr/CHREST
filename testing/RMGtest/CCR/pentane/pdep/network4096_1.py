species(
    label = '[CH]=CC([CH2])[CH][C]=C(16337)',
    structure = SMILES('[CH]=CC([CH2])C=[C][CH2]'),
    E0 = (860.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.6502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155481,'amu*angstrom^2'), symmetry=1, barrier=(3.57481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713687,'amu*angstrom^2'), symmetry=1, barrier=(16.4091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86966,'amu*angstrom^2'), symmetry=1, barrier=(65.9791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797695,0.0724935,-7.75867e-05,4.82855e-08,-1.24602e-11,103629,28.6571], Tmin=(100,'K'), Tmax=(931.974,'K')), NASAPolynomial(coeffs=[10.0712,0.0326914,-1.3525e-05,2.45989e-09,-1.67489e-13,101900,-15.429], Tmin=(931.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C3H3(5450)',
    structure = SMILES('[CH]=C=C'),
    E0 = (338.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,180,180,2603.58],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CHCHCH(4849)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6426,0.0163332,3.86245e-05,-6.71404e-08,2.83615e-11,41729.6,13.2819], Tmin=(100,'K'), Tmax=(937.716,'K')), NASAPolynomial(coeffs=[12.9703,0.0066915,-1.00084e-06,1.67635e-10,-1.71464e-14,38279.7,-43.9468], Tmin=(937.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]C=CC[CH][C]=C(16375)',
    structure = SMILES('[CH]C=CCC=[C][CH2]'),
    E0 = (776.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.625572,0.0641804,-4.1426e-05,1.34585e-08,-1.7774e-12,93548.2,30.4383], Tmin=(100,'K'), Tmax=(1742.09,'K')), NASAPolynomial(coeffs=[15.8261,0.0292786,-1.13743e-05,1.95836e-09,-1.27056e-13,88252.1,-51.3333], Tmin=(1742.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(776.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC[CH]C=[C][CH2](18294)',
    structure = SMILES('[CH]=CC[CH]C=[C][CH2]'),
    E0 = (805.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728647,0.0622572,-4.524e-05,1.66523e-08,-2.46614e-12,96971.9,28.6595], Tmin=(100,'K'), Tmax=(1589.95,'K')), NASAPolynomial(coeffs=[15.7792,0.024393,-9.51785e-06,1.67393e-09,-1.10978e-13,92186,-50.9297], Tmin=(1589.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(805.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(1342)',
    structure = SMILES('C#C'),
    E0 = (218.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,647.284,647.285,3592.16],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08745,0.00578572,8.56332e-06,-1.72824e-08,7.83594e-12,26273.8,4.4608], Tmin=(100,'K'), Tmax=(897.945,'K')), NASAPolynomial(coeffs=[5.89068,0.00209055,4.88462e-08,-5.66903e-11,4.15085e-15,25415.9,-10.7352], Tmin=(897.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.303,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C][CH]C=C(15984)',
    structure = SMILES('[CH2]C=C[C]=C'),
    E0 = (374.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(2.0825,'amu*angstrom^2'), symmetry=1, barrier=(47.8809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07959,'amu*angstrom^2'), symmetry=1, barrier=(47.8139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2251,0.0302874,1.05284e-05,-3.68284e-08,1.72731e-11,45167.7,17.3268], Tmin=(100,'K'), Tmax=(923.512,'K')), NASAPolynomial(coeffs=[10.3879,0.0174193,-5.09535e-06,8.16559e-10,-5.51188e-14,42701,-26.5962], Tmin=(923.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
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
    label = '[CH]=C[CH]C=[C][CH2](17835)',
    structure = SMILES('[CH]C=CC=[C][CH2]'),
    E0 = (748.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,1685,370,261.506,261.508,261.509,261.509],'cm^-1')),
        HinderedRotor(inertia=(1.0436,'amu*angstrom^2'), symmetry=1, barrier=(50.6454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04363,'amu*angstrom^2'), symmetry=1, barrier=(50.6454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04363,'amu*angstrom^2'), symmetry=1, barrier=(50.6454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5943,0.0441971,-6.62019e-06,-2.08649e-08,1.08195e-11,90112.3,23.3821], Tmin=(100,'K'), Tmax=(983.437,'K')), NASAPolynomial(coeffs=[10.7884,0.0265549,-9.84048e-06,1.74264e-09,-1.19723e-13,87348.7,-25.6772], Tmin=(983.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(748.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C(C=[CH])C=[C][CH2](18295)',
    structure = SMILES('[CH]C(C=[CH])C=[C][CH2]'),
    E0 = (1103.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,424.202,424.388,2024.26],'cm^-1')),
        HinderedRotor(inertia=(0.580812,'amu*angstrom^2'), symmetry=1, barrier=(13.354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0178278,'amu*angstrom^2'), symmetry=1, barrier=(2.27817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104466,'amu*angstrom^2'), symmetry=1, barrier=(13.3539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537124,'amu*angstrom^2'), symmetry=1, barrier=(68.639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83956,0.072499,-7.96601e-05,4.8867e-08,-1.23503e-11,132869,27.7325], Tmin=(100,'K'), Tmax=(950.731,'K')), NASAPolynomial(coeffs=[10.7842,0.0306606,-1.36529e-05,2.58361e-09,-1.80272e-13,130978,-19.7427], Tmin=(950.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1103.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=CC([CH2])C=[CH](18296)',
    structure = SMILES('[CH][C]=CC([CH2])C=[CH]'),
    E0 = (1079.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,393.445,396.538,397.803,402.183],'cm^-1')),
        HinderedRotor(inertia=(0.487045,'amu*angstrom^2'), symmetry=1, barrier=(54.3671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.4825,'amu*angstrom^2'), symmetry=1, barrier=(54.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48811,'amu*angstrom^2'), symmetry=1, barrier=(54.2834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486478,'amu*angstrom^2'), symmetry=1, barrier=(54.3495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726533,0.0758135,-8.9224e-05,6.61233e-08,-2.06654e-11,129991,29.737], Tmin=(100,'K'), Tmax=(820.579,'K')), NASAPolynomial(coeffs=[7.77395,0.0387416,-1.64879e-05,2.99269e-09,-2.01849e-13,128926,-2.31147], Tmin=(820.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1079.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=CC([CH2])C=[C][CH2](18297)',
    structure = SMILES('[C]=CC([CH2])C=[C][CH2]'),
    E0 = (1171.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,245.74,3100.57],'cm^-1')),
        HinderedRotor(inertia=(0.237876,'amu*angstrom^2'), symmetry=1, barrier=(10.2054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237765,'amu*angstrom^2'), symmetry=1, barrier=(10.2057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54039,'amu*angstrom^2'), symmetry=1, barrier=(66.1163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54159,'amu*angstrom^2'), symmetry=1, barrier=(66.1167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726207,0.0767157,-0.000103359,8.18817e-08,-2.63436e-11,141034,28.8828], Tmin=(100,'K'), Tmax=(835.356,'K')), NASAPolynomial(coeffs=[8.70915,0.032793,-1.42591e-05,2.61011e-09,-1.76277e-13,139899,-7.00434], Tmin=(835.356,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1171.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=CC1C=C([CH2])C1(18298)',
    structure = SMILES('[CH]=CC1[CH]C(=C)C1'),
    E0 = (576.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7336,0.0336544,4.01949e-05,-7.23671e-08,2.87289e-11,69429.6,22.4154], Tmin=(100,'K'), Tmax=(992.076,'K')), NASAPolynomial(coeffs=[12.9868,0.0278238,-1.0776e-05,2.06096e-09,-1.50945e-13,65250.9,-41.5923], Tmin=(992.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CC1C=CC1(17928)',
    structure = SMILES('[CH2][C]=CC1C=CC1'),
    E0 = (569.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34513,0.046656,4.9418e-07,-3.10765e-08,1.43051e-11,68599,23.0298], Tmin=(100,'K'), Tmax=(1015.98,'K')), NASAPolynomial(coeffs=[12.7982,0.0279923,-1.09685e-05,2.04782e-09,-1.457e-13,64907.8,-39.1185], Tmin=(1015.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC([CH2])C=C1(18299)',
    structure = SMILES('[CH2]C1=CC([CH2])C=C1'),
    E0 = (407.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26125,0.0454361,1.51835e-05,-5.50036e-08,2.53506e-11,49141.2,20.5334], Tmin=(100,'K'), Tmax=(950.784,'K')), NASAPolynomial(coeffs=[15.2513,0.0231615,-7.38898e-06,1.29126e-09,-9.22373e-14,44827.4,-54.9499], Tmin=(950.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CC(=C)C=C[CH2](18300)',
    structure = SMILES('[CH]=CC([CH2])=CC=C'),
    E0 = (473.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763313,0.0558628,-3.77603e-06,-4.42365e-08,2.38587e-11,57050,24.8098], Tmin=(100,'K'), Tmax=(934.474,'K')), NASAPolynomial(coeffs=[18.7687,0.0175713,-4.56069e-06,7.32847e-10,-5.2839e-14,51991.6,-69.8954], Tmin=(934.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CC(=C)C=C(17278)',
    structure = SMILES('[CH2]C(C=C)=C[C]=C'),
    E0 = (425.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.67391,'amu*angstrom^2'), symmetry=1, barrier=(38.4864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67811,'amu*angstrom^2'), symmetry=1, barrier=(38.5831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67848,'amu*angstrom^2'), symmetry=1, barrier=(38.5915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762838,0.0588316,-1.59899e-05,-2.90974e-08,1.81929e-11,51262.2,24.0056], Tmin=(100,'K'), Tmax=(923.691,'K')), NASAPolynomial(coeffs=[17.0045,0.0204652,-5.5988e-06,8.70605e-10,-5.88226e-14,46898,-60.4441], Tmin=(923.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC(=C)C=[C][CH2](18301)',
    structure = SMILES('[CH]=CC([CH2])=C[C]=C'),
    E0 = (672.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.65768,'amu*angstrom^2'), symmetry=1, barrier=(38.1132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66034,'amu*angstrom^2'), symmetry=1, barrier=(38.1746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65863,'amu*angstrom^2'), symmetry=1, barrier=(38.1351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696172,0.0622628,-3.22876e-05,-1.22466e-08,1.25368e-11,80981.6,24.6485], Tmin=(100,'K'), Tmax=(914.848,'K')), NASAPolynomial(coeffs=[17.2418,0.0176387,-4.56898e-06,6.72921e-10,-4.40506e-14,76794.3,-60.0415], Tmin=(914.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C#C[CH2](18302)',
    structure = SMILES('[CH]=CC([CH2])C#C[CH2]'),
    E0 = (763.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2100,2250,500,550,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0477943,'amu*angstrom^2'), symmetry=1, barrier=(11.3502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489413,'amu*angstrom^2'), symmetry=1, barrier=(11.2526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.53129,'amu*angstrom^2'), symmetry=1, barrier=(81.1913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147463,'amu*angstrom^2'), symmetry=1, barrier=(81.0485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896999,0.0673264,-6.78175e-05,3.7977e-08,-8.66059e-12,91896.5,27.5447], Tmin=(100,'K'), Tmax=(1057.78,'K')), NASAPolynomial(coeffs=[11.5813,0.0269226,-1.05208e-05,1.86481e-09,-1.25452e-13,89636.2,-24.6007], Tmin=(1057.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Isobutyl) + radical(Propargyl) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([CH2])C=[C][CH2](18303)',
    structure = SMILES('C#CC([CH2])C=[C][CH2]'),
    E0 = (755.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2175,525,1685,370,750,770,3400,2100,289.313],'cm^-1')),
        HinderedRotor(inertia=(0.183701,'amu*angstrom^2'), symmetry=1, barrier=(10.9112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19245,'amu*angstrom^2'), symmetry=1, barrier=(70.8272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317883,'amu*angstrom^2'), symmetry=1, barrier=(18.8813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19241,'amu*angstrom^2'), symmetry=1, barrier=(70.8271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.652104,0.0738328,-8.32905e-05,5.18258e-08,-1.29743e-11,90941.3,26.105], Tmin=(100,'K'), Tmax=(972.817,'K')), NASAPolynomial(coeffs=[12.1633,0.0265005,-1.03069e-05,1.80974e-09,-1.20647e-13,88701.7,-29.1127], Tmin=(972.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH2](16918)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,260.76,263.05,263.653],'cm^-1')),
        HinderedRotor(inertia=(1.04394,'amu*angstrom^2'), symmetry=1, barrier=(50.8215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2703,0.0129731,8.09839e-06,-1.37802e-08,4.71189e-12,73935.1,11.2848], Tmin=(100,'K'), Tmax=(1104.68,'K')), NASAPolynomial(coeffs=[4.54119,0.0159815,-6.32008e-06,1.15742e-09,-7.99393e-14,73190,2.92522], Tmin=(1104.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C2H2(T)(1343)',
    structure = SMILES('[CH]=[CH]'),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([597.285,1197.51,1197.81,1199.36,2685.83,3758.62],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88059,-0.000843477,2.04095e-05,-2.51946e-08,9.47585e-12,71059.3,4.72114], Tmin=(100,'K'), Tmax=(935.365,'K')), NASAPolynomial(coeffs=[4.9214,0.00358273,-9.24444e-07,1.57274e-10,-1.19541e-14,70476.2,-2.3065], Tmin=(935.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]C=C[CH2](15066)',
    structure = SMILES('[CH]C=C[CH2]'),
    E0 = (493.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,273.892,274.283,274.415],'cm^-1')),
        HinderedRotor(inertia=(0.949581,'amu*angstrom^2'), symmetry=1, barrier=(50.5951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.945904,'amu*angstrom^2'), symmetry=1, barrier=(50.5956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54977,0.0249465,9.89964e-06,-2.57642e-08,1.03662e-11,59460.2,14.975], Tmin=(100,'K'), Tmax=(1021.28,'K')), NASAPolynomial(coeffs=[7.69501,0.020986,-8.0647e-06,1.48637e-09,-1.04578e-13,57564.8,-14.0901], Tmin=(1021.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[CH][CH2](16919)',
    structure = SMILES('[CH2][C]=C[CH][CH2]'),
    E0 = (683.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1114.42],'cm^-1')),
        HinderedRotor(inertia=(0.00740552,'amu*angstrom^2'), symmetry=1, barrier=(6.52457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.14018,'amu*angstrom^2'), symmetry=1, barrier=(95.191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0340784,'amu*angstrom^2'), symmetry=1, barrier=(30.0845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16721,0.0338824,-8.66428e-06,-9.5211e-09,5.10585e-12,82235.7,20.5922], Tmin=(100,'K'), Tmax=(1070.03,'K')), NASAPolynomial(coeffs=[9.10766,0.0208121,-8.38972e-06,1.55219e-09,-1.08387e-13,80013.4,-16.8053], Tmin=(1070.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[C]([CH2])C=[C][CH2](18304)',
    structure = SMILES('[CH]C=C([CH2])C=[C][CH2]'),
    E0 = (862.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,346.9,346.903,346.903,346.905],'cm^-1')),
        HinderedRotor(inertia=(0.594529,'amu*angstrom^2'), symmetry=1, barrier=(50.7709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594516,'amu*angstrom^2'), symmetry=1, barrier=(50.7709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594518,'amu*angstrom^2'), symmetry=1, barrier=(50.7709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594516,'amu*angstrom^2'), symmetry=1, barrier=(50.771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754944,0.0616471,-2.83822e-05,-7.53895e-09,7.6325e-12,103843,27.276], Tmin=(100,'K'), Tmax=(989.655,'K')), NASAPolynomial(coeffs=[13.9382,0.0292088,-1.08118e-05,1.90891e-09,-1.30765e-13,100213,-41.346], Tmin=(989.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(862.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC([CH2])[C]=[C][CH2](18305)',
    structure = SMILES('[CH]=CC([CH2])[C]=[C][CH2]'),
    E0 = (1098.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.15278,'amu*angstrom^2'), symmetry=1, barrier=(3.5127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523546,'amu*angstrom^2'), symmetry=1, barrier=(12.0374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154375,'amu*angstrom^2'), symmetry=1, barrier=(3.54938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7023,'amu*angstrom^2'), symmetry=1, barrier=(62.1312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641333,0.0786428,-0.000108271,8.5967e-08,-2.73954e-11,132237,29.6858], Tmin=(100,'K'), Tmax=(856.682,'K')), NASAPolynomial(coeffs=[9.15267,0.0317156,-1.3522e-05,2.44151e-09,-1.6316e-13,131043,-8.52081], Tmin=(856.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1098.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH2])C=[C][CH2](18306)',
    structure = SMILES('[CH]=[C]C([CH2])C=[C][CH2]'),
    E0 = (1098.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.15278,'amu*angstrom^2'), symmetry=1, barrier=(3.5127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523546,'amu*angstrom^2'), symmetry=1, barrier=(12.0374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154375,'amu*angstrom^2'), symmetry=1, barrier=(3.54938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7023,'amu*angstrom^2'), symmetry=1, barrier=(62.1312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641333,0.0786428,-0.000108271,8.5967e-08,-2.73954e-11,132237,29.6858], Tmin=(100,'K'), Tmax=(856.682,'K')), NASAPolynomial(coeffs=[9.15267,0.0317156,-1.3522e-05,2.44151e-09,-1.6316e-13,131043,-8.52081], Tmin=(856.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1098.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C](C)C=[C][CH2](18307)',
    structure = SMILES('[CH]C=C(C)C=[C][CH2]'),
    E0 = (710.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656089,0.0652558,-4.0554e-05,1.00685e-08,2.25029e-14,85624.2,27.3535], Tmin=(100,'K'), Tmax=(1123.23,'K')), NASAPolynomial(coeffs=[12.5176,0.0338258,-1.3018e-05,2.29381e-09,-1.54176e-13,82277.6,-34.2856], Tmin=(1123.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(710.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC([CH2])[C]=C[CH2](18308)',
    structure = SMILES('[CH]=CC([CH2])[C]=C[CH2]'),
    E0 = (860.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.6502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155481,'amu*angstrom^2'), symmetry=1, barrier=(3.57481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713687,'amu*angstrom^2'), symmetry=1, barrier=(16.4091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86966,'amu*angstrom^2'), symmetry=1, barrier=(65.9791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797695,0.0724935,-7.75867e-05,4.82855e-08,-1.24602e-11,103629,28.6571], Tmin=(100,'K'), Tmax=(931.974,'K')), NASAPolynomial(coeffs=[10.0712,0.0326914,-1.3525e-05,2.45989e-09,-1.67489e-13,101900,-15.429], Tmin=(931.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CC([CH2])[C]=C(17281)',
    structure = SMILES('[CH2][C]=CC([CH2])[C]=C'),
    E0 = (851.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1670,1700,300,440,367.224,367.423],'cm^-1')),
        HinderedRotor(inertia=(0.0746498,'amu*angstrom^2'), symmetry=1, barrier=(7.14554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746324,'amu*angstrom^2'), symmetry=1, barrier=(7.14599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746136,'amu*angstrom^2'), symmetry=1, barrier=(7.14575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.653593,'amu*angstrom^2'), symmetry=1, barrier=(62.5797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01168,0.0710235,-7.38545e-05,3.90021e-08,-4.86818e-12,102505,27.9886], Tmin=(100,'K'), Tmax=(649.865,'K')), NASAPolynomial(coeffs=[8.45328,0.0353862,-1.50639e-05,2.76471e-09,-1.88633e-13,101323,-6.35725], Tmin=(649.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(851.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)[C]=[C][CH2](18309)',
    structure = SMILES('[CH]=CC(C)[C]=[C][CH2]'),
    E0 = (893.433,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,231.6],'cm^-1')),
        HinderedRotor(inertia=(0.249762,'amu*angstrom^2'), symmetry=1, barrier=(9.57599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249369,'amu*angstrom^2'), symmetry=1, barrier=(9.57407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250913,'amu*angstrom^2'), symmetry=1, barrier=(9.57391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.758684,'amu*angstrom^2'), symmetry=1, barrier=(28.9269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827832,0.0743868,-8.6051e-05,5.9882e-08,-1.75802e-11,107565,27.2459], Tmin=(100,'K'), Tmax=(817.206,'K')), NASAPolynomial(coeffs=[8.57472,0.0364716,-1.64636e-05,3.11899e-09,-2.16985e-13,106299,-8.56535], Tmin=(817.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(893.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(C)C=[C][CH2](18310)',
    structure = SMILES('[CH]=[C]C(C)C=[C][CH2]'),
    E0 = (893.433,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,231.6],'cm^-1')),
        HinderedRotor(inertia=(0.249762,'amu*angstrom^2'), symmetry=1, barrier=(9.57599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249369,'amu*angstrom^2'), symmetry=1, barrier=(9.57407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250913,'amu*angstrom^2'), symmetry=1, barrier=(9.57391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.758684,'amu*angstrom^2'), symmetry=1, barrier=(28.9269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827832,0.0743868,-8.6051e-05,5.9882e-08,-1.75802e-11,107565,27.2459], Tmin=(100,'K'), Tmax=(817.206,'K')), NASAPolynomial(coeffs=[8.57472,0.0364716,-1.64636e-05,3.11899e-09,-2.16985e-13,106299,-8.56535], Tmin=(817.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(893.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]([CH2])C=C[CH2](18311)',
    structure = SMILES('[CH]C=C([CH2])C=C[CH2]'),
    E0 = (624.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801231,0.0568893,-3.03529e-06,-3.75416e-08,1.89146e-11,75239.5,26.6362], Tmin=(100,'K'), Tmax=(966.254,'K')), NASAPolynomial(coeffs=[15.0442,0.0298497,-1.06151e-05,1.87897e-09,-1.31015e-13,70996.8,-49.3003], Tmin=(966.254,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[C]([CH2])C=C(17280)',
    structure = SMILES('[CH2][C]=CC([CH2])=C[CH2]'),
    E0 = (609.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,562.923],'cm^-1')),
        HinderedRotor(inertia=(0.284772,'amu*angstrom^2'), symmetry=1, barrier=(63.8987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284826,'amu*angstrom^2'), symmetry=1, barrier=(63.8368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284774,'amu*angstrom^2'), symmetry=1, barrier=(63.8838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283281,'amu*angstrom^2'), symmetry=1, barrier=(63.843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05238,0.0526186,-4.8284e-06,-3.37637e-08,1.77357e-11,73451,26.6575], Tmin=(100,'K'), Tmax=(952.27,'K')), NASAPolynomial(coeffs=[14.6495,0.0251847,-8.36739e-06,1.44442e-09,-1.00284e-13,69515.6,-45.3419], Tmin=(952.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])[C]=[C]C(18312)',
    structure = SMILES('[CH]=CC([CH2])[C]=[C]C'),
    E0 = (947.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,184.015],'cm^-1')),
        HinderedRotor(inertia=(0.329438,'amu*angstrom^2'), symmetry=1, barrier=(7.91699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329377,'amu*angstrom^2'), symmetry=1, barrier=(7.91703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329492,'amu*angstrom^2'), symmetry=1, barrier=(7.91676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329355,'amu*angstrom^2'), symmetry=1, barrier=(7.91733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.636342,0.081223,-0.000117274,1.00201e-07,-3.3941e-11,114014,29.422], Tmin=(100,'K'), Tmax=(864.737,'K')), NASAPolynomial(coeffs=[6.88667,0.0377221,-1.65101e-05,3.00777e-09,-2.01408e-13,113478,3.33006], Tmin=(864.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(947.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH2])C=C[CH2](18313)',
    structure = SMILES('[CH]=[C]C([CH2])C=C[CH2]'),
    E0 = (860.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.6502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155481,'amu*angstrom^2'), symmetry=1, barrier=(3.57481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.713687,'amu*angstrom^2'), symmetry=1, barrier=(16.4091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.86966,'amu*angstrom^2'), symmetry=1, barrier=(65.9791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797695,0.0724935,-7.75867e-05,4.82855e-08,-1.24602e-11,103629,28.6571], Tmin=(100,'K'), Tmax=(931.974,'K')), NASAPolynomial(coeffs=[10.0712,0.0326914,-1.3525e-05,2.45989e-09,-1.67489e-13,101900,-15.429], Tmin=(931.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])C=C(17282)',
    structure = SMILES('[CH2][C]=[C]C([CH2])C=C'),
    E0 = (851.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1670,1700,300,440,367.224,367.423],'cm^-1')),
        HinderedRotor(inertia=(0.0746498,'amu*angstrom^2'), symmetry=1, barrier=(7.14554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746324,'amu*angstrom^2'), symmetry=1, barrier=(7.14599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0746136,'amu*angstrom^2'), symmetry=1, barrier=(7.14575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.653593,'amu*angstrom^2'), symmetry=1, barrier=(62.5797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01168,0.0710235,-7.38545e-05,3.90021e-08,-4.86818e-12,102505,27.9886], Tmin=(100,'K'), Tmax=(649.865,'K')), NASAPolynomial(coeffs=[8.45328,0.0353862,-1.50639e-05,2.76471e-09,-1.88633e-13,101323,-6.35725], Tmin=(649.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(851.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[C]([CH2])C=[C]C(18314)',
    structure = SMILES('[CH]C=C([CH2])C=[C]C'),
    E0 = (744.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347661,0.0719253,-5.58085e-05,2.28158e-08,-3.80199e-12,89657.6,27.1873], Tmin=(100,'K'), Tmax=(1417.1,'K')), NASAPolynomial(coeffs=[15.1391,0.0301742,-1.16152e-05,2.02543e-09,-1.34239e-13,85465.4,-49.3295], Tmin=(1417.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C([CH2])C=[C]C(18315)',
    structure = SMILES('[CH]=[C]C([CH2])C=[C]C'),
    E0 = (947.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,184.015],'cm^-1')),
        HinderedRotor(inertia=(0.329438,'amu*angstrom^2'), symmetry=1, barrier=(7.91699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329377,'amu*angstrom^2'), symmetry=1, barrier=(7.91703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329492,'amu*angstrom^2'), symmetry=1, barrier=(7.91676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329355,'amu*angstrom^2'), symmetry=1, barrier=(7.91733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.636342,0.081223,-0.000117274,1.00201e-07,-3.3941e-11,114014,29.422], Tmin=(100,'K'), Tmax=(864.737,'K')), NASAPolynomial(coeffs=[6.88667,0.0377221,-1.65101e-05,3.00777e-09,-2.01408e-13,113478,3.33006], Tmin=(864.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(947.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([C]=C)[CH][CH2](18316)',
    structure = SMILES('[CH]=CC([C]=C)[CH][CH2]'),
    E0 = (916.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,233.105,2831.38],'cm^-1')),
        HinderedRotor(inertia=(0.00304614,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158414,'amu*angstrom^2'), symmetry=1, barrier=(6.84602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175202,'amu*angstrom^2'), symmetry=1, barrier=(6.81632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.736104,'amu*angstrom^2'), symmetry=1, barrier=(28.8727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63456,0.0618315,-3.40097e-05,-4.79913e-08,6.3895e-11,110260,28.5228], Tmin=(100,'K'), Tmax=(484.489,'K')), NASAPolynomial(coeffs=[6.37388,0.0397799,-1.86072e-05,3.56549e-09,-2.48754e-13,109600,7.02428], Tmin=(484.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(916.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1CC1[C]=C(18317)',
    structure = SMILES('[CH]=CC1CC1[C]=C'),
    E0 = (671.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4237,0.0398573,2.88128e-05,-6.84223e-08,2.96966e-11,80872,25.7791], Tmin=(100,'K'), Tmax=(961.246,'K')), NASAPolynomial(coeffs=[15.8043,0.021438,-7.08267e-06,1.30235e-09,-9.6581e-14,76193.6,-52.9847], Tmin=(961.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1C=CC1[C]=C(18318)',
    structure = SMILES('[CH2]C1C=CC1[C]=C'),
    E0 = (627.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36304,0.0474315,-3.12439e-06,-2.87347e-08,1.43716e-11,75556.1,23.9651], Tmin=(100,'K'), Tmax=(974.286,'K')), NASAPolynomial(coeffs=[12.3942,0.0268644,-9.52099e-06,1.68629e-09,-1.17187e-13,72233.3,-34.9878], Tmin=(974.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(=C)C[C]=C(16332)',
    structure = SMILES('[CH]=CC(=C)C[C]=C'),
    E0 = (623.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.775977,'amu*angstrom^2'), symmetry=1, barrier=(17.8412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773011,'amu*angstrom^2'), symmetry=1, barrier=(17.773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.774472,'amu*angstrom^2'), symmetry=1, barrier=(17.8066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627552,0.0651991,-4.57238e-05,9.70892e-09,1.83306e-12,75166.4,27.2354], Tmin=(100,'K'), Tmax=(1020.59,'K')), NASAPolynomial(coeffs=[15.6221,0.0231853,-8.59921e-06,1.54391e-09,-1.07109e-13,71233.2,-49.6851], Tmin=(1020.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(C)=C[C]=C(18319)',
    structure = SMILES('[CH]C=C(C)C=C=C'),
    E0 = (531.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434752,0.0693313,-4.83658e-05,1.53364e-08,-1.31543e-12,64068.8,25.1719], Tmin=(100,'K'), Tmax=(1151.6,'K')), NASAPolynomial(coeffs=[14.3857,0.0312023,-1.21551e-05,2.16223e-09,-1.46244e-13,60170.8,-47.0765], Tmin=(1151.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC([CH2])[CH]C#C(18320)',
    structure = SMILES('[CH]=CC([CH2])[CH]C#C'),
    E0 = (796.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2175,525,3120,650,792.5,1650,750,770,3400,2100,215.639],'cm^-1')),
        HinderedRotor(inertia=(1.87498,'amu*angstrom^2'), symmetry=1, barrier=(58.1689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00361267,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00350044,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79423,'amu*angstrom^2'), symmetry=1, barrier=(58.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.523311,0.0687477,-6.98162e-05,3.85562e-08,-8.28574e-12,95921.4,29.3783], Tmin=(100,'K'), Tmax=(1263.19,'K')), NASAPolynomial(coeffs=[13.9111,0.0203171,-5.13756e-06,6.37644e-10,-3.24372e-14,93020.8,-36.4315], Tmin=(1263.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(796.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC=C[C]=C(17856)',
    structure = SMILES('[CH]C=CC=C=C'),
    E0 = (569.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05172,'amu*angstrom^2'), symmetry=1, barrier=(47.1731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05189,'amu*angstrom^2'), symmetry=1, barrier=(47.177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38512,0.0481531,-1.41309e-05,-1.5813e-08,9.49606e-12,68556.4,21.1552], Tmin=(100,'K'), Tmax=(995.318,'K')), NASAPolynomial(coeffs=[12.4634,0.0242384,-9.14602e-06,1.64941e-09,-1.14882e-13,65330.4,-37.3668], Tmin=(995.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[C]([CH2])C[C]=C(16336)',
    structure = SMILES('[CH]C=C([CH2])C[C]=C'),
    E0 = (773.695,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.585982,0.0669476,-4.64228e-05,1.65959e-08,-2.43261e-12,93183.5,29.5164], Tmin=(100,'K'), Tmax=(1576.56,'K')), NASAPolynomial(coeffs=[14.7357,0.0310476,-1.22664e-05,2.15255e-09,-1.42302e-13,88721.9,-45.1897], Tmin=(1576.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(773.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]C([CH2])C=[CH](14729)',
    structure = SMILES('[CH]C=CC([CH2])C=[CH]'),
    E0 = (842.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650,381.893,381.931,381.96,382.054],'cm^-1')),
        HinderedRotor(inertia=(0.518375,'amu*angstrom^2'), symmetry=1, barrier=(53.6676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.51829,'amu*angstrom^2'), symmetry=1, barrier=(53.6649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518315,'amu*angstrom^2'), symmetry=1, barrier=(53.6627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518223,'amu*angstrom^2'), symmetry=1, barrier=(53.6681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.866593,0.0699301,-5.97814e-05,3.04148e-08,-6.69257e-12,101383,28.7616], Tmin=(100,'K'), Tmax=(1058.24,'K')), NASAPolynomial(coeffs=[9.01011,0.0391488,-1.61506e-05,2.92844e-09,-1.99181e-13,99659.7,-10.9872], Tmin=(1058.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(842.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C([CH2])C[C]=C(16338)',
    structure = SMILES('[CH]=[C]C([CH2])C[C]=C'),
    E0 = (959.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,373.179,373.575],'cm^-1')),
        HinderedRotor(inertia=(0.00218815,'amu*angstrom^2'), symmetry=1, barrier=(11.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11921,'amu*angstrom^2'), symmetry=1, barrier=(11.8354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118622,'amu*angstrom^2'), symmetry=1, barrier=(11.8247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.761934,'amu*angstrom^2'), symmetry=1, barrier=(75.6195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938518,0.0680787,-6.68991e-05,3.78605e-08,-8.91034e-12,115548,31.1661], Tmin=(100,'K'), Tmax=(1015.51,'K')), NASAPolynomial(coeffs=[10.248,0.0314095,-1.27352e-05,2.30264e-09,-1.56609e-13,113657,-13.8898], Tmin=(1015.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(959.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC([CH2])C=[CH](16339)',
    structure = SMILES('[CH]=[C]CC([CH2])C=[CH]'),
    E0 = (969.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3115,3125,620,680,785,800,1600,1700,282.837],'cm^-1')),
        HinderedRotor(inertia=(0.132955,'amu*angstrom^2'), symmetry=1, barrier=(7.57683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133508,'amu*angstrom^2'), symmetry=1, barrier=(7.57566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210597,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489359,'amu*angstrom^2'), symmetry=1, barrier=(27.8656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778217,0.0684668,-6.41125e-05,3.31395e-08,-6.99132e-12,116670,31.6665], Tmin=(100,'K'), Tmax=(1137.47,'K')), NASAPolynomial(coeffs=[12.2659,0.0280694,-1.08396e-05,1.91637e-09,-1.28891e-13,114056,-25.2345], Tmin=(1137.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(969.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(C)C=[CH](18321)',
    structure = SMILES('[CH][C]=CC(C)C=[CH]'),
    E0 = (874.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,285.09,285.539,287.333,287.946],'cm^-1')),
        HinderedRotor(inertia=(0.89535,'amu*angstrom^2'), symmetry=1, barrier=(52.5691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.876628,'amu*angstrom^2'), symmetry=1, barrier=(52.5586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886902,'amu*angstrom^2'), symmetry=1, barrier=(52.5576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895448,'amu*angstrom^2'), symmetry=1, barrier=(52.5637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.970591,0.070872,-6.45534e-05,3.67142e-08,-9.32765e-12,105317,27.0909], Tmin=(100,'K'), Tmax=(910.77,'K')), NASAPolynomial(coeffs=[7.21293,0.0434559,-1.93996e-05,3.662e-09,-2.54925e-13,104180,-2.44119], Tmin=(910.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(874.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CC([CH2])C=C(17267)',
    structure = SMILES('[CH][C]=CC([CH2])C=C'),
    E0 = (832.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.923291,0.0707368,-6.66183e-05,4.02263e-08,-1.07258e-11,100266,28.6345], Tmin=(100,'K'), Tmax=(879.702,'K')), NASAPolynomial(coeffs=[7.21221,0.0421403,-1.78563e-05,3.27186e-09,-2.23474e-13,99159.9,-0.89965], Tmin=(879.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(832.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    E0 = (860.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1020.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1020.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (860.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1186.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1315.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1291.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1383.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (868.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (868.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (867.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (924.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (924.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (901.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (990.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (987.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (974.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (978.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (860.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (928.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1108.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1273.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1074.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1310.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1310.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (959.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1064.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (966.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1045.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1045.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (1054.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1052.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1108.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1254.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (904.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (991.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (980.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1079.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (866.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (868.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (883.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (883.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (1023.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (990.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (1018.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (1055.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (1101.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (1114.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (907.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (1051.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['C3H3(5450)', 'CH2CHCHCH(4849)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]C=CC[CH][C]=C(16375)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=CC[CH]C=[C][CH2](18294)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['C2H2(1342)', 'C=[C][CH]C=C(15984)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH]=C[CH]C=[C][CH2](17835)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C(C=[CH])C=[C][CH2](18295)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH][C]=CC([CH2])C=[CH](18296)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[C]=CC([CH2])C=[C][CH2](18297)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=CC1C=C([CH2])C1(18298)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2][C]=CC1C=CC1(17928)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2]C1=CC([CH2])C=C1(18299)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_DSSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=CC(=C)C=C[CH2](18300)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2][C]=CC(=C)C=C(17278)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]=CC(=C)C=[C][CH2](18301)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-TwoDe_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH]=CC([CH2])C#C[CH2](18302)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C#CC([CH2])C=[C][CH2](18303)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.276e+10,'cm^3/(mol*s)'), n=1.103, Ea=(20.5476,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 139 used for Ct-Cs_Ct-H;HJ
Exact match found for rate rule [Ct-Cs_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C][CH2](16918)', 'CH2CHCHCH(4849)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00509465,'m^3/(mol*s)'), n=2.41843, Ea=(13.2728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C2H2(T)(1343)', 'C=[C][CH]C=C(15984)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0101893,'m^3/(mol*s)'), n=2.41843, Ea=(13.2728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C3H3(5450)', '[CH]C=C[CH2](15066)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(28.2068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 22.9 to 28.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H2(1342)', '[CH2][C]=C[CH][CH2](16919)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C][CH2](16918)', '[CH]C=C[CH2](15066)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C2H2(T)(1343)', '[CH2][C]=C[CH][CH2](16919)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH]=C[C]([CH2])C=[C][CH2](18304)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH]=CC([CH2])[C]=[C][CH2](18305)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]=[C]C([CH2])C=[C][CH2](18306)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=C[C](C)C=[C][CH2](18307)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC([CH2])[C]=C[CH2](18308)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2][C]=CC([CH2])[C]=C(17281)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=CC(C)[C]=[C][CH2](18309)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]C(C)C=[C][CH2](18310)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=C[C]([CH2])C=C[CH2](18311)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2][C]=C[C]([CH2])C=C(17280)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC([CH2])[C]=[C]C(18312)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]C([CH2])C=C[CH2](18313)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;XH_out] for rate rule [R4H_SSD;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2][C]=[C]C([CH2])C=C(17282)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=C[C]([CH2])C=[C]C(18314)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]C([CH2])C=[C]C(18315)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC([C]=C)[CH][CH2](18316)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=CC1CC1[C]=C(18317)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH2]C1C=CC1[C]=C(18318)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=CC(=C)C[C]=C(16332)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=CC(C)=C[C]=C(18319)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)', '[CH]=CC([CH2])[CH]C#C(18320)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CH2(T)(20)', '[CH]=CC=C[C]=C(17856)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=C[C]([CH2])C[C]=C(16336)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH]=C[CH]C([CH2])C=[CH](14729)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C]C([CH2])C[C]=C(16338)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=[C]CC([CH2])C=[CH](16339)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=[C][CH]C(C)C=[CH](18321)'],
    products = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=CC([CH2])[CH][C]=C(16337)'],
    products = ['[CH][C]=CC([CH2])C=C(17267)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6Hall;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4096',
    isomers = [
        '[CH]=CC([CH2])[CH][C]=C(16337)',
    ],
    reactants = [
        ('C3H3(5450)', 'CH2CHCHCH(4849)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4096',
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


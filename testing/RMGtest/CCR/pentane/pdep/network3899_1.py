species(
    label = 'C=[C][CH]C(=C)[O](15778)',
    structure = SMILES('[CH2][C]=CC(=C)[O]'),
    E0 = (337.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,331.383,335.37],'cm^-1')),
        HinderedRotor(inertia=(0.251608,'amu*angstrom^2'), symmetry=1, barrier=(19.8263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.940183,'amu*angstrom^2'), symmetry=1, barrier=(73.9665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16324,0.0521249,-4.82109e-05,2.30793e-08,-4.26991e-12,40727.4,22.7626], Tmin=(100,'K'), Tmax=(1461.18,'K')), NASAPolynomial(coeffs=[13.7022,0.0131394,-3.40593e-06,4.54312e-10,-2.54632e-14,37560.5,-40.7837], Tmin=(1461.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[CH2][C]=C[C]=C(17191)',
    structure = SMILES('[CH2][C]=C[C]=C'),
    E0 = (612.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(2.17267,'amu*angstrom^2'), symmetry=1, barrier=(49.9539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.16922,'amu*angstrom^2'), symmetry=1, barrier=(49.8747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20155,0.0347757,-1.38759e-05,-8.03191e-09,6.4918e-12,73770.6,17.1921], Tmin=(100,'K'), Tmax=(915.397,'K')), NASAPolynomial(coeffs=[9.16602,0.0169746,-5.40492e-06,8.73114e-10,-5.70733e-14,71966.3,-18.6825], Tmin=(915.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C)"""),
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
    label = '[CH][C]=CC(=C)[O](17192)',
    structure = SMILES('[CH]=[C]C=C([CH2])[O]'),
    E0 = (586.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.04542,'amu*angstrom^2'), symmetry=1, barrier=(24.0363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04507,'amu*angstrom^2'), symmetry=1, barrier=(24.0283,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908408,0.0615872,-7.58456e-05,4.57077e-08,-1.03049e-11,70614.5,23.1374], Tmin=(100,'K'), Tmax=(1237.86,'K')), NASAPolynomial(coeffs=[15.1626,0.00746192,-4.85966e-07,-1.41344e-10,1.7735e-14,67703.4,-46.1766], Tmin=(1237.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=CC(=C)O1(17193)',
    structure = SMILES('C=C1[CH]C(=C)O1'),
    E0 = (188.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32341,0.0174386,6.49104e-05,-1.00128e-07,4.00727e-11,22766.6,16.9515], Tmin=(100,'K'), Tmax=(959.586,'K')), NASAPolynomial(coeffs=[15.252,0.0124858,-3.84855e-06,7.90756e-10,-6.62219e-14,18032.2,-56.6285], Tmin=(959.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2][C]=[C]C([CH2])[O](17194)',
    structure = SMILES('[CH2][C]=[C]C([CH2])[O]'),
    E0 = (843.779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,180,809.451],'cm^-1')),
        HinderedRotor(inertia=(0.0681146,'amu*angstrom^2'), symmetry=1, barrier=(1.56609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809771,'amu*angstrom^2'), symmetry=1, barrier=(18.6182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0400344,'amu*angstrom^2'), symmetry=1, barrier=(18.6205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44711,0.0559746,-5.90488e-05,3.25705e-08,-7.21201e-12,101575,25.9871], Tmin=(100,'K'), Tmax=(1091.34,'K')), NASAPolynomial(coeffs=[11.5221,0.0190475,-8.29385e-06,1.56578e-09,-1.09545e-13,99376.2,-23.4995], Tmin=(1091.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(843.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]1CO1(17195)',
    structure = SMILES('C=[C][CH][C]1CO1'),
    E0 = (511.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66523,0.0501254,-3.81719e-05,3.8854e-09,7.52414e-12,61623.4,18.9074], Tmin=(100,'K'), Tmax=(760.169,'K')), NASAPolynomial(coeffs=[10.7202,0.0177049,-4.2438e-06,4.80237e-10,-2.1729e-14,59806.8,-25.1881], Tmin=(760.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]1[CH]C(=O)C1(16429)',
    structure = SMILES('[CH2][C]1C=C([O])C1'),
    E0 = (373.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38279,0.0200451,5.00902e-05,-8.65424e-08,3.73737e-11,44937.2,19.3924], Tmin=(100,'K'), Tmax=(908.176,'K')), NASAPolynomial(coeffs=[13.867,0.0101016,-6.06221e-07,-5.68593e-11,3.12843e-15,41175.4,-44.1325], Tmin=(908.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Allyl_T) + radical(Isobutyl)"""),
)

species(
    label = '[O][C]1C=[C]CC1(17196)',
    structure = SMILES('[O]C1=C[C]CC1'),
    E0 = (393.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26808,0.021057,5.02221e-05,-8.56119e-08,3.5781e-11,47400.1,18.3449], Tmin=(100,'K'), Tmax=(941.045,'K')), NASAPolynomial(coeffs=[14.8452,0.0108096,-2.32404e-06,4.10408e-10,-3.5333e-14,43119.6,-51.7346], Tmin=(941.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC1([CH2])[O](17197)',
    structure = SMILES('[CH2]C1([O])[CH]C1=C'),
    E0 = (524.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50969,0.0442824,-1.31213e-05,-1.92917e-08,1.17031e-11,63215.6,19.5051], Tmin=(100,'K'), Tmax=(975.125,'K')), NASAPolynomial(coeffs=[14.6645,0.0147278,-5.20262e-06,9.62506e-10,-7.03414e-14,59489.7,-49.5778], Tmin=(975.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(C=CCJCO) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1([O])C=[C]C1(17198)',
    structure = SMILES('[CH2]C1([O])C=[C]C1'),
    E0 = (617.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03951,0.0340801,3.17012e-07,-2.45168e-08,1.16586e-11,74391.5,22.0024], Tmin=(100,'K'), Tmax=(995.741,'K')), NASAPolynomial(coeffs=[11.2345,0.0179771,-6.81045e-06,1.26821e-09,-9.09799e-14,71527.4,-27.5054], Tmin=(995.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C#CC(=C)[O](17199)',
    structure = SMILES('[CH2]C#CC(=C)[O]'),
    E0 = (308.845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,2100,2250,500,550,372.799,373.083],'cm^-1')),
        HinderedRotor(inertia=(0.290481,'amu*angstrom^2'), symmetry=1, barrier=(28.7111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290355,'amu*angstrom^2'), symmetry=1, barrier=(28.7027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73037,0.0436397,-3.08645e-05,5.33772e-09,2.05795e-12,37232.8,20.2564], Tmin=(100,'K'), Tmax=(996.028,'K')), NASAPolynomial(coeffs=[12.3892,0.0136304,-4.9415e-06,8.85011e-10,-6.18442e-14,34474.7,-34.3099], Tmin=(996.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=C(C)OJ) + radical(Propargyl)"""),
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
    label = '[CH2][C]=[C]C(=C)[O](17200)',
    structure = SMILES('[CH2][C]=[C]C(=C)[O]'),
    E0 = (536.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,201.63,204.896],'cm^-1')),
        HinderedRotor(inertia=(2.07625,'amu*angstrom^2'), symmetry=1, barrier=(61.2376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10266,'amu*angstrom^2'), symmetry=1, barrier=(61.2859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38752,0.0552231,-6.57979e-05,4.17797e-08,-1.03038e-11,64646.2,21.5467], Tmin=(100,'K'), Tmax=(1080.9,'K')), NASAPolynomial(coeffs=[11.3537,0.0145903,-4.20398e-06,5.79245e-10,-3.18889e-14,62710.8,-26.2958], Tmin=(1080.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])[CH][C]=C(16203)',
    structure = SMILES('[CH]C([O])=C[C]=C'),
    E0 = (550.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.16134,'amu*angstrom^2'), symmetry=1, barrier=(49.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15849,'amu*angstrom^2'), symmetry=1, barrier=(49.6278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.973876,0.0578595,-5.6507e-05,2.89348e-08,-5.76231e-12,66364.9,22.9252], Tmin=(100,'K'), Tmax=(1341.27,'K')), NASAPolynomial(coeffs=[13.7717,0.0154074,-4.23793e-06,5.72638e-10,-3.18357e-14,63317.4,-41.1376], Tmin=(1341.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=[C]C(=C)[O](17201)',
    structure = SMILES('[CH2]C=[C]C(=C)[O]'),
    E0 = (298.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14801,0.0538039,-5.19828e-05,2.65345e-08,-5.23314e-12,36054.7,21.9344], Tmin=(100,'K'), Tmax=(1384.84,'K')), NASAPolynomial(coeffs=[13.0084,0.0142678,-3.44147e-06,4.14172e-10,-2.08513e-14,33275.9,-37.3191], Tmin=(1384.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=[C]C(=C)O(17202)',
    structure = SMILES('[CH2][C]=[C]C(=C)O'),
    E0 = (398.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(1.61576,'amu*angstrom^2'), symmetry=1, barrier=(37.1494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62627,'amu*angstrom^2'), symmetry=1, barrier=(37.3912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61916,'amu*angstrom^2'), symmetry=1, barrier=(37.2277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778862,0.0620733,-6.84183e-05,3.79407e-08,-7.96884e-12,48099.6,22.2285], Tmin=(100,'K'), Tmax=(1313.75,'K')), NASAPolynomial(coeffs=[15.3301,0.0107772,-1.86718e-06,1.18227e-10,-5.38201e-16,44879.6,-49.6474], Tmin=(1313.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)[CH][C]=C(16205)',
    structure = SMILES('[CH]C(O)=C[C]=C'),
    E0 = (413.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05963,'amu*angstrom^2'), symmetry=1, barrier=(47.3549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06202,'amu*angstrom^2'), symmetry=1, barrier=(47.41,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0611,'amu*angstrom^2'), symmetry=1, barrier=(47.3887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.206949,0.0664197,-6.44406e-05,3.11387e-08,-5.67404e-12,49825.6,24.1861], Tmin=(100,'K'), Tmax=(1551.76,'K')), NASAPolynomial(coeffs=[17.3329,0.0121232,-2.14336e-06,1.59157e-10,-3.87597e-15,45732.6,-62.0241], Tmin=(1551.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([O])[C]=[C]C(17203)',
    structure = SMILES('C=C([O])[C]=[C]C'),
    E0 = (418.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.800719,'amu*angstrom^2'), symmetry=1, barrier=(18.4101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800325,'amu*angstrom^2'), symmetry=1, barrier=(18.401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4401,0.0596822,-7.05754e-05,3.82572e-08,-4.4622e-12,50440.5,19.8314], Tmin=(100,'K'), Tmax=(680.289,'K')), NASAPolynomial(coeffs=[10.2952,0.0192506,-7.08132e-06,1.17629e-09,-7.46543e-14,48966.5,-21.4571], Tmin=(680.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C([O])C=C[CH2](9868)',
    structure = SMILES('[CH]=C([O])C=C[CH2]'),
    E0 = (346.959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.01783,'amu*angstrom^2'), symmetry=1, barrier=(23.4019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.69404,'amu*angstrom^2'), symmetry=1, barrier=(61.9412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5367,0.0462158,-2.33927e-05,-1.02198e-08,9.849e-12,41825.7,21.3484], Tmin=(100,'K'), Tmax=(914.094,'K')), NASAPolynomial(coeffs=[14.0776,0.0126023,-3.12849e-06,4.50296e-10,-2.94381e-14,38644.6,-42.8871], Tmin=(914.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C=[C]C(17204)',
    structure = SMILES('[CH]=C([O])C=[C]C'),
    E0 = (466.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.678072,'amu*angstrom^2'), symmetry=1, barrier=(15.5902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675272,'amu*angstrom^2'), symmetry=1, barrier=(15.5258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25697,0.0593345,-7.01754e-05,4.35589e-08,-1.06379e-11,56236.2,21.2681], Tmin=(100,'K'), Tmax=(1006.02,'K')), NASAPolynomial(coeffs=[12.0816,0.0162923,-5.99437e-06,1.02491e-09,-6.74069e-14,54058.3,-31.0187], Tmin=(1006.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH2][C]=C[C]=O(9897)',
    structure = SMILES('[CH2][C]=C[C]=O'),
    E0 = (432.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.00399,'amu*angstrom^2'), symmetry=1, barrier=(23.0837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00442,'amu*angstrom^2'), symmetry=1, barrier=(23.0936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33986,0.0429161,-4.87447e-05,1.23757e-08,1.41321e-11,52082.9,16.5479], Tmin=(100,'K'), Tmax=(512.282,'K')), NASAPolynomial(coeffs=[6.6704,0.019476,-1.04855e-05,2.11536e-09,-1.51263e-13,51503.1,-2.77645], Tmin=(512.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'H2CC(T)(1341)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.3669e-08,8.45456e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.859e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C([CH2])[O](9653)',
    structure = SMILES('[CH]C(=C)[O]'),
    E0 = (322.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,387.097,387.511,388.12,388.534],'cm^-1')),
        HinderedRotor(inertia=(0.483854,'amu*angstrom^2'), symmetry=1, barrier=(51.3097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72438,0.0238014,-4.51958e-06,-8.49326e-09,4.40325e-12,38798.3,14.6668], Tmin=(100,'K'), Tmax=(1019.83,'K')), NASAPolynomial(coeffs=[7.15237,0.0155204,-5.70458e-06,1.01811e-09,-7.00768e-14,37422.6,-9.09941], Tmin=(1019.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C=C1CO1(17205)',
    structure = SMILES('C=[C]C=C1CO1'),
    E0 = (275.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60692,0.0370652,1.5896e-05,-6.13122e-08,3.07079e-11,33181.9,16.2955], Tmin=(100,'K'), Tmax=(906.893,'K')), NASAPolynomial(coeffs=[18.4869,0.00410674,1.77898e-06,-4.83739e-10,3.19724e-14,28413.9,-72.8985], Tmin=(906.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1C=C([O])C1(16199)',
    structure = SMILES('C=C1C=C([O])C1'),
    E0 = (142.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00816,0.0242072,5.14884e-05,-9.63447e-08,4.23359e-11,17285.3,17.2213], Tmin=(100,'K'), Tmax=(920.7,'K')), NASAPolynomial(coeffs=[18.3514,0.00397859,1.72241e-06,-4.11476e-10,2.26176e-14,12123.8,-71.9624], Tmin=(920.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=C=C=C(17206)',
    structure = SMILES('[CH2]C#CC(=C)O'),
    E0 = (171.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38691,0.0473427,-2.23832e-05,-1.30398e-08,1.06571e-11,20674.7,19.9883], Tmin=(100,'K'), Tmax=(951.253,'K')), NASAPolynomial(coeffs=[15.9848,0.0105116,-3.02228e-06,5.25363e-10,-3.90604e-14,16786.6,-55.5478], Tmin=(951.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Propargyl)"""),
)

species(
    label = 'C=C=C=C(C)[O](17207)',
    structure = SMILES('[CH2]C#CC(C)=O'),
    E0 = (153.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40325,0.0362264,-2.07599e-05,5.76663e-09,-6.59962e-13,18501.3,19.3597], Tmin=(100,'K'), Tmax=(1859.82,'K')), NASAPolynomial(coeffs=[8.94468,0.0221576,-9.41309e-06,1.69931e-09,-1.13232e-13,16068.1,-16.2578], Tmin=(1859.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cs-CtHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C]1OC1[C]=C(17208)',
    structure = SMILES('[CH2][C]1OC1[C]=C'),
    E0 = (593.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19906,0.0532682,-5.49809e-05,3.05568e-08,-6.42895e-12,71471.3,23.8565], Tmin=(100,'K'), Tmax=(1370.67,'K')), NASAPolynomial(coeffs=[11.5419,0.0140607,-2.19804e-06,8.08449e-11,5.73753e-15,69483.7,-26.2103], Tmin=(1370.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C[C]1[O](17209)',
    structure = SMILES('C=[C]C1C[C]1[O]'),
    E0 = (586.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84973,0.0384765,-9.34036e-06,-1.75656e-08,1.00451e-11,70646.9,21.4127], Tmin=(100,'K'), Tmax=(981.046,'K')), NASAPolynomial(coeffs=[12.4101,0.0157379,-5.64064e-06,1.03191e-09,-7.39832e-14,67597.1,-34.3162], Tmin=(981.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C#CC=C([CH2])[O](17210)',
    structure = SMILES('[CH]=C=CC(=C)[O]'),
    E0 = (312.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.31068,'amu*angstrom^2'), symmetry=1, barrier=(30.1352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987814,0.0569827,-6.40073e-05,3.53993e-08,-7.32832e-12,37748.2,21.0311], Tmin=(100,'K'), Tmax=(1355.48,'K')), NASAPolynomial(coeffs=[15.0497,0.00757548,-5.78155e-07,-1.02623e-10,1.36391e-14,34662.8,-48.406], Tmin=(1355.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C=C=C[C]=O(9908)',
    structure = SMILES('C=C=C[C]=O'),
    E0 = (213.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.13978,'amu*angstrom^2'), symmetry=1, barrier=(26.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00684,0.0284132,-1.9291e-05,6.01948e-09,-7.90367e-13,25757,13.5223], Tmin=(100,'K'), Tmax=(1551.73,'K')), NASAPolynomial(coeffs=[6.76337,0.0187298,-9.93048e-06,1.99793e-09,-1.42457e-13,24591.2,-6.25135], Tmin=(1551.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = 'C=[C]C=C=C(17211)',
    structure = SMILES('C=[C]C=C=C'),
    E0 = (433.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.5987,'amu*angstrom^2'), symmetry=1, barrier=(36.7571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0932,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98554,0.0388285,-2.1811e-05,-2.2868e-09,4.80017e-12,52215,14.9887], Tmin=(100,'K'), Tmax=(940.235,'K')), NASAPolynomial(coeffs=[10.8095,0.0147077,-4.73738e-06,7.85956e-10,-5.27182e-14,49962.5,-30.1922], Tmin=(940.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[CH]C(=C)[O](14143)',
    structure = SMILES('[CH]C=CC(=C)[O]'),
    E0 = (352.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15066,'amu*angstrom^2'), symmetry=1, barrier=(49.4479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.14942,'amu*angstrom^2'), symmetry=1, barrier=(49.4195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32006,0.0516581,-3.01759e-05,-1.32923e-09,5.53682e-12,42498.1,21.2726], Tmin=(100,'K'), Tmax=(951.664,'K')), NASAPolynomial(coeffs=[13.0056,0.0196593,-6.72024e-06,1.13999e-09,-7.69582e-14,39498.8,-38.597], Tmin=(951.664,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C][C]=C(C)[O](17212)',
    structure = SMILES('C=[C][C]=C(C)[O]'),
    E0 = (379.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.79342,'amu*angstrom^2'), symmetry=1, barrier=(41.2342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00796,'amu*angstrom^2'), symmetry=1, barrier=(23.1751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12706,0.0617989,-7.68488e-05,5.09037e-08,-1.30609e-11,45703.8,20.861], Tmin=(100,'K'), Tmax=(1050.99,'K')), NASAPolynomial(coeffs=[11.7842,0.0162579,-4.74325e-06,6.56622e-10,-3.60284e-14,43738.7,-29.775], Tmin=(1050.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]C=C([CH2])O(17213)',
    structure = SMILES('[CH]=[C]C=C([CH2])O'),
    E0 = (448.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.12158,'amu*angstrom^2'), symmetry=1, barrier=(25.7874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11934,'amu*angstrom^2'), symmetry=1, barrier=(25.7357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12287,'amu*angstrom^2'), symmetry=1, barrier=(25.8171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.235894,0.0691565,-8.08321e-05,4.47381e-08,-9.10843e-12,54070.8,24.0507], Tmin=(100,'K'), Tmax=(1408.54,'K')), NASAPolynomial(coeffs=[18.6594,0.00434233,1.49639e-06,-5.25933e-10,4.31763e-14,50120.2,-66.7434], Tmin=(1408.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C=C(C)[O](17214)',
    structure = SMILES('[CH]=[C]C=C(C)[O]'),
    E0 = (427.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.896911,'amu*angstrom^2'), symmetry=1, barrier=(20.6217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896065,'amu*angstrom^2'), symmetry=1, barrier=(20.6023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26488,0.0571076,-5.8098e-05,2.64811e-08,-3.0413e-12,51485.7,21.1792], Tmin=(100,'K'), Tmax=(873.884,'K')), NASAPolynomial(coeffs=[13.2881,0.0138265,-3.98e-06,5.852e-10,-3.56346e-14,48935.6,-37.7728], Tmin=(873.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C1OC1=C(17174)',
    structure = SMILES('C=[C]C1OC1=C'),
    E0 = (320.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74711,0.0326184,2.51541e-05,-6.66155e-08,3.06702e-11,38590.4,17.6265], Tmin=(100,'K'), Tmax=(939.329,'K')), NASAPolynomial(coeffs=[18.1591,0.00563033,-2.55744e-07,3.94856e-11,-1.02113e-14,33614.6,-70.5991], Tmin=(939.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])C[C]=C(15562)',
    structure = SMILES('[CH]=C([O])C[C]=C'),
    E0 = (496.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,339.344,340.743],'cm^-1')),
        HinderedRotor(inertia=(0.113631,'amu*angstrom^2'), symmetry=1, barrier=(9.29722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113113,'amu*angstrom^2'), symmetry=1, barrier=(9.30356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.78,'J/mol'), sigma=(6.06207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.87 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41567,0.0518712,-5.0883e-05,2.59077e-08,-5.20182e-12,59864.6,24.535], Tmin=(100,'K'), Tmax=(1216.07,'K')), NASAPolynomial(coeffs=[12.6633,0.0148747,-5.24865e-06,8.9038e-10,-5.87625e-14,57129,-31.9286], Tmin=(1216.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(=C)[O](15779)',
    structure = SMILES('[CH]=[C]CC(=C)[O]'),
    E0 = (496.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,339.344,340.743],'cm^-1')),
        HinderedRotor(inertia=(0.113631,'amu*angstrom^2'), symmetry=1, barrier=(9.29722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113113,'amu*angstrom^2'), symmetry=1, barrier=(9.30356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41567,0.0518712,-5.0883e-05,2.59077e-08,-5.20182e-12,59864.6,24.535], Tmin=(100,'K'), Tmax=(1216.07,'K')), NASAPolynomial(coeffs=[12.6633,0.0148747,-5.24865e-06,8.9038e-10,-5.87625e-14,57129,-31.9286], Tmin=(1216.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH]C[C]=O(15745)',
    structure = SMILES('[CH2][C]=CC[C]=O'),
    E0 = (420.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950,1130.06],'cm^-1')),
        HinderedRotor(inertia=(0.0482167,'amu*angstrom^2'), symmetry=1, barrier=(17.6152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76571,'amu*angstrom^2'), symmetry=1, barrier=(17.6052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76638,'amu*angstrom^2'), symmetry=1, barrier=(17.6206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.61,'J/mol'), sigma=(5.8495,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.13 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78239,0.04201,-2.54007e-05,5.06591e-09,1.61991e-13,50616.2,23.0246], Tmin=(100,'K'), Tmax=(1352.83,'K')), NASAPolynomial(coeffs=[12.9962,0.0178357,-8.55618e-06,1.67302e-09,-1.17996e-13,46760.2,-37.5023], Tmin=(1352.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]1C=C([CH2])O1(17215)',
    structure = SMILES('[CH2][C]1[CH]C(=C)O1'),
    E0 = (435.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0599847,0.0654185,-6.75562e-05,3.30814e-08,-5.90046e-12,52603.3,20.5459], Tmin=(100,'K'), Tmax=(1662.63,'K')), NASAPolynomial(coeffs=[18.1137,0.00532763,1.42383e-06,-4.98642e-10,3.90892e-14,48822.4,-69.5678], Tmin=(1662.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C2CsJOC(O)) + radical(C=CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1C=[C]CO1(17216)',
    structure = SMILES('[CH2]C1=C[C]CO1'),
    E0 = (436.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14161,0.0242225,4.61647e-05,-8.95797e-08,4.06483e-11,52619.5,15.6117], Tmin=(100,'K'), Tmax=(893.007,'K')), NASAPolynomial(coeffs=[16.2698,0.0056744,2.17699e-06,-6.43444e-10,4.56731e-14,48312.5,-60.9374], Tmin=(893.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(2,3-Dihydrofuran) + radical(C=C(O)CJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=[C]C1CC1=O(17170)',
    structure = SMILES('C=[C]C1CC1=O'),
    E0 = (312.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80736,0.0430217,-3.18392e-05,1.24575e-08,-1.99222e-12,37724.3,19.0449], Tmin=(100,'K'), Tmax=(1467.89,'K')), NASAPolynomial(coeffs=[10.6747,0.018858,-7.14676e-06,1.24296e-09,-8.22347e-14,35121.1,-27.1386], Tmin=(1467.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CC(=C)[O](15777)',
    structure = SMILES('C=C=CC(=C)[O]'),
    E0 = (158.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13355,'amu*angstrom^2'), symmetry=1, barrier=(26.0625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976076,0.0558197,-5.47945e-05,2.69e-08,-5.05269e-12,19170.5,20.4568], Tmin=(100,'K'), Tmax=(1432.29,'K')), NASAPolynomial(coeffs=[15.4408,0.010735,-2.6683e-06,3.52107e-10,-1.99502e-14,15507.9,-52.8453], Tmin=(1432.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C][O](2821)',
    structure = SMILES('[CH2][C][O]'),
    E0 = (648.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,2091.34],'cm^-1')),
        HinderedRotor(inertia=(0.0328816,'amu*angstrom^2'), symmetry=1, barrier=(10.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26468,0.0187769,-3.12134e-05,3.2382e-08,-1.33361e-11,78037.6,11.2589], Tmin=(100,'K'), Tmax=(748.797,'K')), NASAPolynomial(coeffs=[3.86755,0.0107985,-5.69995e-06,1.18131e-09,-8.60475e-14,78080.7,9.41548], Tmin=(748.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C=[C][C]([CH2])[O](17217)',
    structure = SMILES('[CH2][CH][C]=C([CH2])[O]'),
    E0 = (613.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,2902.9],'cm^-1')),
        HinderedRotor(inertia=(0.0126662,'amu*angstrom^2'), symmetry=1, barrier=(15.3129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0656804,'amu*angstrom^2'), symmetry=1, barrier=(79.3704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0584521,'amu*angstrom^2'), symmetry=1, barrier=(15.3064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45754,0.0521996,-5.26766e-05,2.82778e-08,-6.02221e-12,73919.8,25.1033], Tmin=(100,'K'), Tmax=(1146.31,'K')), NASAPolynomial(coeffs=[11.7488,0.0162882,-5.68405e-06,9.47564e-10,-6.16223e-14,71560.4,-25.9509], Tmin=(1146.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=C(O)CJ) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH][C]([CH2])[O](17218)',
    structure = SMILES('[CH][CH]C=C([CH2])[O]'),
    E0 = (618.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,180.357,827.147,827.456,828.328],'cm^-1')),
        HinderedRotor(inertia=(0.0306035,'amu*angstrom^2'), symmetry=1, barrier=(14.8061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132221,'amu*angstrom^2'), symmetry=1, barrier=(78.3634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00407491,'amu*angstrom^2'), symmetry=1, barrier=(1.97142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45519,0.0492725,-3.46515e-05,3.42501e-09,4.29078e-12,74539.4,23.1633], Tmin=(100,'K'), Tmax=(936.103,'K')), NASAPolynomial(coeffs=[13.7056,0.0135693,-4.10995e-06,6.66918e-10,-4.49173e-14,71516.7,-39.0238], Tmin=(936.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=C(O)CJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])[C]=C=C(17219)',
    structure = SMILES('[CH2]C#CC([CH2])[O]'),
    E0 = (532.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2100,2250,500,550,946.759,946.814],'cm^-1')),
        HinderedRotor(inertia=(0.972415,'amu*angstrom^2'), symmetry=1, barrier=(22.3577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259321,'amu*angstrom^2'), symmetry=1, barrier=(16.4956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.34679,'amu*angstrom^2'), symmetry=1, barrier=(76.9493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57004,0.0466647,-2.98403e-05,6.24882e-10,4.65977e-12,64174.4,24.6021], Tmin=(100,'K'), Tmax=(950.958,'K')), NASAPolynomial(coeffs=[12.9817,0.0147187,-4.77422e-06,8.05775e-10,-5.50368e-14,61278,-33.6957], Tmin=(950.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CC(C)OJ) + radical(CJCO) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=CC([CH2])[O](17220)',
    structure = SMILES('[CH]=C=CC([CH2])[O]'),
    E0 = (547.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,540,610,2055,3120,650,792.5,1650,188.034],'cm^-1')),
        HinderedRotor(inertia=(0.825185,'amu*angstrom^2'), symmetry=1, barrier=(20.4659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82321,'amu*angstrom^2'), symmetry=1, barrier=(20.4682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0806,0.0556284,-5.49764e-05,2.75105e-08,-5.32939e-12,65982.7,25.7538], Tmin=(100,'K'), Tmax=(1334.73,'K')), NASAPolynomial(coeffs=[14.789,0.0122228,-3.5851e-06,5.37562e-10,-3.2963e-14,62530.3,-43.5637], Tmin=(1334.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CJCO) + radical(C=C=CJ)"""),
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
    E0 = (337.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1132.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (797.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (345.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (866.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (524.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (475.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (454.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (524.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (617.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (535.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (599.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (522.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (775.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (748.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (762.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (487.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (549.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (605.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (611.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (750.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (499.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (870.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (977.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (343.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (345.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (415.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (415.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (593.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (586.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (534.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (634.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (691.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (520.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (535.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (481.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (460.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (340.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (689.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (541.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (665.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (475.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (436.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (343.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (337.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (1042.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (636.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (641.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (710.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (862.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['ketene(1375)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH2][C]=C[C]=C(17191)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.15242e+13,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH][C]=CC(=C)[O](17192)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2]C1=CC(=C)O1(17193)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][C]=[C]C([CH2])[O](17194)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2][C]=C[C]1CO1(17195)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2][C]1[CH]C(=O)C1(16429)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[O][C]1C=[C]CC1(17196)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.49459e+10,'s^-1'), n=0.314867, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2]C1=CC1([CH2])[O](17197)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.69117e+12,'s^-1'), n=0.1675, Ea=(187.075,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra_2H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 186.9 to 187.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2]C1([O])C=[C]C1(17198)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.36716e+08,'s^-1'), n=1.01412, Ea=(280.168,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 278.2 to 280.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]C#CC(=C)[O](17199)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['ketene(1375)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=O(1376)', 'C3H3(5450)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=O(1376)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][C]=[C]C(=C)[O](17200)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]=C([O])[CH][C]=C(16203)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2]C=[C]C(=C)[O](17201)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.44274e+08,'s^-1'), n=1.26608, Ea=(149.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=[C]C(=C)O(17202)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.87605e+06,'s^-1'), n=1.79171, Ea=(150.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)[CH][C]=C(16205)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])[C]=[C]C(17203)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C([O])C=C[CH2](9868)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C([O])C=[C]C(17204)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(T)(20)', '[CH2][C]=C[C]=O(9897)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H2CC(T)(1341)', '[CH]=C([CH2])[O](9653)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=[C]C=C1CO1(17205)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=C1C=C([O])C1(16199)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2]C(O)=C=C=C(17206)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=C=C=C(C)[O](17207)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2][C]1OC1[C]=C(17208)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(255.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=[C]C1C[C]1[O](17209)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(248.977,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 247.7 to 249.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', 'C#CC=C([CH2])[O](17210)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.156e+09,'cm^3/(mol*s)'), n=1.502, Ea=(9.92026,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 194 used for Ct-H_Ct-Cd;HJ
Exact match found for rate rule [Ct-H_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2(T)(20)', 'C=C=C[C]=O(9908)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(4)', 'C=[C]C=C=C(17211)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.40891,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;O_atom_triplet]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C[CH]C(=C)[O](14143)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C][C]=C(C)[O](17212)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=[C]C=C([CH2])O(17213)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]C=C(C)[O](17214)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=[C]C1OC1=C(17174)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C([O])C[C]=C(15562)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC(=C)[O](15779)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C][CH]C[C]=O(15745)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2][C]1C=C([CH2])O1(17215)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['[CH2][C]1C=[C]CO1(17216)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.19156e+09,'s^-1'), n=0.640131, Ea=(99.1115,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 96.3 to 99.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=[C]C1CC1=O(17170)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C][CH]C(=C)[O](15778)'],
    products = ['C=C=CC(=C)[O](15777)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][C][O](2821)', 'C3H3(5450)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C=[C][C]([CH2])[O](17217)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C[CH][C]([CH2])[O](17218)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([O])[C]=C=C(17219)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=C=CC([CH2])[O](17220)'],
    products = ['C=[C][CH]C(=C)[O](15778)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3899',
    isomers = [
        'C=[C][CH]C(=C)[O](15778)',
    ],
    reactants = [
        ('ketene(1375)', 'C3H3(5450)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3899',
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


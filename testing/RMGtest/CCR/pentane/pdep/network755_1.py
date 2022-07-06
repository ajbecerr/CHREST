species(
    label = 'C#COC(=C)O(2352)',
    structure = SMILES('C#COC(=C)O'),
    E0 = (22.0036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05372,'amu*angstrom^2'), symmetry=1, barrier=(24.2271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05374,'amu*angstrom^2'), symmetry=1, barrier=(24.2275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05325,'amu*angstrom^2'), symmetry=1, barrier=(24.2162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16789,0.0604999,-7.40004e-05,4.47468e-08,-1.04797e-11,2750,18.4496], Tmin=(100,'K'), Tmax=(1053.49,'K')), NASAPolynomial(coeffs=[14.0065,0.0117524,-4.59094e-06,8.22806e-10,-5.61352e-14,44.9662,-44.158], Tmin=(1053.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.0036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH)"""),
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
    label = 'O[C]1CC=[C]O1(3073)',
    structure = SMILES('O[C]1CC=[C]O1'),
    E0 = (116.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45484,0.0190253,4.46186e-05,-8.14961e-08,3.62259e-11,14049.6,18.8629], Tmin=(100,'K'), Tmax=(899.12,'K')), NASAPolynomial(coeffs=[14.5354,0.00502379,1.67472e-06,-4.93634e-10,3.39284e-14,10270.8,-47.0677], Tmin=(899.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = 'C#COC([CH2])[O](2845)',
    structure = SMILES('C#COC([CH2])[O]'),
    E0 = (260.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,339.582,339.679,339.842],'cm^-1')),
        HinderedRotor(inertia=(0.298483,'amu*angstrom^2'), symmetry=1, barrier=(24.4313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298424,'amu*angstrom^2'), symmetry=1, barrier=(24.4313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298468,'amu*angstrom^2'), symmetry=1, barrier=(24.4355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21649,0.0590493,-6.65808e-05,3.67522e-08,-7.914e-12,31479.1,20.7875], Tmin=(100,'K'), Tmax=(1137.03,'K')), NASAPolynomial(coeffs=[14.3101,0.0129863,-5.81255e-06,1.12205e-09,-7.98686e-14,28501.5,-44.063], Tmin=(1137.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = '[C]#COC([CH2])O(3074)',
    structure = SMILES('[C]#COC([CH2])O'),
    E0 = (372.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2175,525,329.414,329.414,329.414],'cm^-1')),
        HinderedRotor(inertia=(0.258952,'amu*angstrom^2'), symmetry=1, barrier=(19.9402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258952,'amu*angstrom^2'), symmetry=1, barrier=(19.9402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258952,'amu*angstrom^2'), symmetry=1, barrier=(19.9402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258952,'amu*angstrom^2'), symmetry=1, barrier=(19.9402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796232,0.0663341,-8.05938e-05,4.38867e-08,-8.15412e-12,44899.1,22.1754], Tmin=(100,'K'), Tmax=(917.004,'K')), NASAPolynomial(coeffs=[17.221,0.00632062,-1.45334e-06,1.84074e-10,-1.08338e-14,41397.7,-58.3085], Tmin=(917.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCO) + radical(Acetyl)"""),
)

species(
    label = '[C]#CO[C](C)O(3075)',
    structure = SMILES('[C]#CO[C](C)O'),
    E0 = (365.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2175,525,364.127,364.376,365.882],'cm^-1')),
        HinderedRotor(inertia=(0.194503,'amu*angstrom^2'), symmetry=1, barrier=(18.5032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195867,'amu*angstrom^2'), symmetry=1, barrier=(18.4879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195839,'amu*angstrom^2'), symmetry=1, barrier=(18.4753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194564,'amu*angstrom^2'), symmetry=1, barrier=(18.4702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975894,0.0620033,-7.17749e-05,3.76679e-08,-6.76081e-12,44130.1,21.5841], Tmin=(100,'K'), Tmax=(940.876,'K')), NASAPolynomial(coeffs=[16.2277,0.00775965,-2.19126e-06,3.34668e-10,-2.17974e-14,40791,-53.5603], Tmin=(940.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cs_P) + radical(Acetyl)"""),
)

species(
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2COH(1378)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989113,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16243,0.0134242,5.5647e-06,-1.95528e-08,9.36447e-12,12455.2,10.1543], Tmin=(100,'K'), Tmax=(925.606,'K')), NASAPolynomial(coeffs=[8.19867,0.00453476,-8.93529e-07,1.26102e-10,-9.46676e-15,10971.4,-16.7326], Tmin=(925.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HC2(2881)',
    structure = SMILES('[C]#C'),
    E0 = (556.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0139907,-3.08144e-05,3.10836e-08,-1.10946e-11,66983,5.75946], Tmin=(100,'K'), Tmax=(918.723,'K')), NASAPolynomial(coeffs=[3.14385,0.00498487,-2.32624e-06,4.08757e-10,-2.55787e-14,67315.6,7.08554], Tmin=(918.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.809,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""HC2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2COOH(2279)',
    structure = SMILES('[CH2]C(=O)O'),
    E0 = (-247.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,412.68,412.705,413.038,413.69],'cm^-1')),
        HinderedRotor(inertia=(0.000988091,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000988507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92688,0.0179295,4.49428e-06,-1.67263e-08,6.86763e-12,-29754.3,12.035], Tmin=(100,'K'), Tmax=(1058,'K')), NASAPolynomial(coeffs=[7.98363,0.0118142,-5.27086e-06,1.04332e-09,-7.61676e-14,-31552.1,-16.0853], Tmin=(1058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2COOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'OH(5)',
    structure = SMILES('[OH]'),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133396,-4.70039e-06,5.64374e-09,-2.06315e-12,3411.96,1.99787], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35656e-07,1.40239e-11,6.34502e-16,3669.56,5.59056], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C#CO[C]=C(3076)',
    structure = SMILES('C#CO[C]=C'),
    E0 = (418.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21813,'amu*angstrom^2'), symmetry=1, barrier=(28.0073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21862,'amu*angstrom^2'), symmetry=1, barrier=(28.0184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30607,0.0360715,-3.58533e-05,1.8608e-08,-3.87254e-12,50442.2,17.7059], Tmin=(100,'K'), Tmax=(1159.65,'K')), NASAPolynomial(coeffs=[9.19036,0.0123252,-5.13738e-06,9.49691e-10,-6.5685e-14,48845.5,-16.5265], Tmin=(1159.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CJO)"""),
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
    label = 'C#COC([CH2])=O(2691)',
    structure = SMILES('C#COC([CH2])=O'),
    E0 = (90.4817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180,180,180,756.309,756.522,756.53],'cm^-1')),
        HinderedRotor(inertia=(0.0964546,'amu*angstrom^2'), symmetry=1, barrier=(39.1512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964858,'amu*angstrom^2'), symmetry=1, barrier=(39.1537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964372,'amu*angstrom^2'), symmetry=1, barrier=(39.155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78591,0.0433021,-3.72489e-05,1.55609e-08,-2.566e-12,10966.6,19.6326], Tmin=(100,'K'), Tmax=(1454.22,'K')), NASAPolynomial(coeffs=[12.9294,0.0126506,-5.63242e-06,1.06676e-09,-7.42587e-14,7725.54,-38.3014], Tmin=(1454.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.4817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCO)"""),
)

species(
    label = '[CH]=C(O)OC#C(3077)',
    structure = SMILES('[CH]=C(O)OC#C'),
    E0 = (269.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3120,650,792.5,1650,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.08182,'amu*angstrom^2'), symmetry=1, barrier=(24.8731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07923,'amu*angstrom^2'), symmetry=1, barrier=(24.8135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07932,'amu*angstrom^2'), symmetry=1, barrier=(24.8157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10947,0.0638236,-8.98893e-05,6.10556e-08,-1.5915e-11,32469.1,19.0638], Tmin=(100,'K'), Tmax=(951.383,'K')), NASAPolynomial(coeffs=[14.1132,0.00915205,-3.69341e-06,6.56714e-10,-4.40057e-14,29994.8,-43.0239], Tmin=(951.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_P)"""),
)

species(
    label = '[C]#COC(=C)O(3078)',
    structure = SMILES('[C]#COC(=C)O'),
    E0 = (359.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,2175,525,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.58119,'amu*angstrom^2'), symmetry=1, barrier=(36.3547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.737426,'amu*angstrom^2'), symmetry=1, barrier=(16.9549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.736636,'amu*angstrom^2'), symmetry=1, barrier=(16.9367,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30326,0.062452,-8.95143e-05,5.96747e-08,-1.35182e-11,43289.8,18.7216], Tmin=(100,'K'), Tmax=(715.443,'K')), NASAPolynomial(coeffs=[12.3008,0.0114034,-4.36984e-06,7.2709e-10,-4.56007e-14,41449,-32.5193], Tmin=(715.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = 'ethynol(3079)',
    structure = SMILES('C#CO'),
    E0 = (79.2211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.11152,0.0160255,-1.46236e-05,6.90897e-09,-1.25524e-12,9563.09,7.82738], Tmin=(100,'K'), Tmax=(1511.92,'K')), NASAPolynomial(coeffs=[7.00147,0.0039905,-9.53776e-07,1.18632e-10,-6.31686e-15,8586.11,-11.8884], Tmin=(1511.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.2211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ethynol""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]C(O)OC#C(3080)',
    structure = SMILES('[CH]C(O)OC#C'),
    E0 = (274.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2175,525,750,770,3400,2100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.843607,0.0648384,-7.80922e-05,4.43801e-08,-9.63696e-12,33155,20.771], Tmin=(100,'K'), Tmax=(1140.91,'K')), NASAPolynomial(coeffs=[17.0999,0.00784387,-3.15866e-06,5.94017e-10,-4.23264e-14,29445.6,-59.7991], Tmin=(1140.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsJ2_singlet-CsH) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = '[C]=COC(=C)O(3081)',
    structure = SMILES('[C]=COC(=C)O'),
    E0 = (198.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.572941,'amu*angstrom^2'), symmetry=1, barrier=(13.173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.573211,'amu*angstrom^2'), symmetry=1, barrier=(13.1793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5733,'amu*angstrom^2'), symmetry=1, barrier=(13.1813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03102,0.0725102,-0.000124657,1.09001e-07,-3.63944e-11,23990.1,20.0323], Tmin=(100,'K'), Tmax=(872.881,'K')), NASAPolynomial(coeffs=[9.13834,0.0205228,-9.82511e-06,1.82677e-09,-1.22242e-13,23140,-14.7412], Tmin=(872.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdOsH) + group(CdJ2_singlet-Cds)"""),
)

species(
    label = 'C#CCC(=O)O(3082)',
    structure = SMILES('C#CCC(=O)O'),
    E0 = (-198.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21421,0.0300301,3.88313e-06,-2.44805e-08,1.04305e-11,-23855.5,18.2559], Tmin=(100,'K'), Tmax=(1059.46,'K')), NASAPolynomial(coeffs=[11.075,0.0174831,-7.95294e-06,1.59332e-09,-1.17198e-13,-26906.5,-30.5419], Tmin=(1059.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CtHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#COC(C)=O(2353)',
    structure = SMILES('C#COC(C)=O'),
    E0 = (-121.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93679,0.0381596,-1.78175e-05,-2.94976e-09,3.21485e-12,-14485.5,17.7624], Tmin=(100,'K'), Tmax=(1119.54,'K')), NASAPolynomial(coeffs=[11.3365,0.0176341,-7.81319e-06,1.5116e-09,-1.07946e-13,-17408.6,-32.3023], Tmin=(1119.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Ct-CtOs) + group(Ct-CtH)"""),
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
    E0 = (142.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (172.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (283.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (397.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (400.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (268.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (309.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (447.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (302.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (480.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (570.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (172.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (331.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (232.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (335.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (226.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#COC(=C)O(2352)'],
    products = ['ketene(1375)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.29914e+17,'s^-1'), n=-1.73308, Ea=(120.08,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.16692063000437535, var=8.814681701569807, Tref=1000.0, N=65, correlation='Root',), comment="""BM rule fitted to 2 training reactions at node Root
    Total Standard Deviation in ln(k): 6.371362717619001
Exact match found for rate rule [Root]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O[C]1CC=[C]O1(3073)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#COC([CH2])[O](2845)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]#COC([CH2])O(3074)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]#CO[C](C)O(3075)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HCCO(2227)', 'CH2COH(1378)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HC2(2881)', 'CH2COOH(2279)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(5)', 'C#CO[C]=C(3076)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=4.95181e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_2R->C_Ext-2C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_2R->C_Ext-2C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_2R->C_Ext-2C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C#COC([CH2])=O(2691)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]=C(O)OC#C(3077)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[C]#COC(=C)O(3078)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.81e+08,'m^3/(mol*s)'), n=-8.9089e-10, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_N-Sp-3C-2CN',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_N-Sp-3C-2CN
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_N-Sp-3C-2CN]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#COC(=C)O(2352)'],
    products = ['ethynol(3079)', 'ketene(1375)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(371062,'s^-1'), n=1.69565, Ea=(150.411,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, correlation='Root_N-1R!H->C_N-3R!H-inRing_N-5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R!H->C_N-3R!H-inRing_N-5R!H->O
    Total Standard Deviation in ln(k): 11.540182761524994
Exact match found for rate rule [Root_N-1R!H->C_N-3R!H-inRing_N-5R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Retroene"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(O)OC#C(3080)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH]
Euclidian distance = 1.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[C]=COC(=C)O(3081)'],
    products = ['C#COC(=C)O(2352)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.99854e+12,'s^-1'), n=-0.201832, Ea=(34.0789,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH] for rate rule [CdJ2=C;CdJ2;CH]
Euclidian distance = 1.4142135623730951
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#COC(=C)O(2352)'],
    products = ['C#CCC(=O)O(3082)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#COC(=C)O(2352)'],
    products = ['C#COC(C)=O(2353)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = 'PDepNetwork #755',
    isomers = [
        'C#COC(=C)O(2352)',
    ],
    reactants = [
        ('ketene(1375)', 'ketene(1375)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #755',
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


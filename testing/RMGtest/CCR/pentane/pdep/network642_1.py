species(
    label = 'vinoxyperoxy(1385)',
    structure = SMILES('[O]OCC=O'),
    E0 = (-95.4746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.546537,'amu*angstrom^2'), symmetry=1, barrier=(12.566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.546366,'amu*angstrom^2'), symmetry=1, barrier=(12.562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3712.02,'J/mol'), sigma=(5.93388,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=579.81 K, Pc=40.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29854,0.0327818,-2.77436e-05,1.17559e-08,-1.97285e-12,-11417.8,17.2937], Tmin=(100,'K'), Tmax=(1433.5,'K')), NASAPolynomial(coeffs=[10.4206,0.010118,-4.02818e-06,7.26608e-10,-4.93432e-14,-13746.4,-24.8153], Tmin=(1433.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.4746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""vinoxyperoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'HO2(9)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26737e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.914,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11273e-07,6.34955e-11,-4.86383e-15,83.4204,3.0934], Tmin=(923.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36761e-06,-4.93131e-09,1.45956e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.7,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00338e-07,1.59031e-10,-1.14892e-14,-1048.44,6.08305], Tmin=(1087.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = 'CO(2039)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CH3OO(44)',
    structure = SMILES('CO[O]'),
    E0 = (1.92253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1236.37,3898.56],'cm^-1')),
        HinderedRotor(inertia=(1.07614,'amu*angstrom^2'), symmetry=1, barrier=(24.7426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0333,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.13,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63773,0.00381152,2.1454e-05,-2.68366e-08,9.49389e-12,247.919,9.81744], Tmin=(100,'K'), Tmax=(1001.47,'K')), NASAPolynomial(coeffs=[4.74153,0.0100949,-3.97179e-06,7.49789e-10,-5.38568e-14,-509.344,1.81365], Tmin=(1001.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.92253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3OO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'OCH2CHO(1527)',
    structure = SMILES('[O]CC=O'),
    E0 = (-90.1149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,584.193],'cm^-1')),
        HinderedRotor(inertia=(0.069063,'amu*angstrom^2'), symmetry=1, barrier=(16.7258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04365,0.0136818,1.69887e-05,-3.22643e-08,1.32999e-11,-10797.3,13.3581], Tmin=(100,'K'), Tmax=(973.599,'K')), NASAPolynomial(coeffs=[8.68695,0.00860874,-3.1005e-06,5.99601e-10,-4.54045e-14,-12754.6,-18.1251], Tmin=(973.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.1149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]C[CH]O[O](2040)',
    structure = SMILES('[O]C[CH]O[O]'),
    E0 = (221.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,180,1849.86],'cm^-1')),
        HinderedRotor(inertia=(0.18576,'amu*angstrom^2'), symmetry=1, barrier=(4.27099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185775,'amu*angstrom^2'), symmetry=1, barrier=(4.27134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27381,0.0497176,-0.000101209,1.06992e-07,-4.0818e-11,26679.1,19.2469], Tmin=(100,'K'), Tmax=(875.886,'K')), NASAPolynomial(coeffs=[0.61055,0.0281083,-1.4187e-05,2.68808e-09,-1.81444e-13,28090.7,33.4459], Tmin=(875.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH]1COOO1(2041)',
    structure = SMILES('[CH]1COOO1'),
    E0 = (138.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72563,0.0185668,2.76988e-05,-5.97715e-08,2.90323e-11,16665.3,12.6763], Tmin=(100,'K'), Tmax=(846.404,'K')), NASAPolynomial(coeffs=[11.5519,0.00485035,2.39344e-06,-7.61707e-10,6.03955e-14,14168.4,-34.3574], Tmin=(846.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(123trioxolane) + radical(CCsJOO)"""),
)

species(
    label = '[O]C1COO1(2042)',
    structure = SMILES('[O]C1COO1'),
    E0 = (9.76301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63211,0.0191775,-4.55024e-06,-1.60941e-09,5.22299e-13,1177.3,10.4054], Tmin=(100,'K'), Tmax=(2084.65,'K')), NASAPolynomial(coeffs=[12.4936,0.0118129,-6.18643e-06,1.13176e-09,-7.24136e-14,-4611.68,-43.879], Tmin=(2084.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.76301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(12dioxetane) + radical(CCOJ)"""),
)

species(
    label = 'HCO(1372)',
    structure = SMILES('[CH]=O'),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[O](46)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (205.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2028.66,2028.67],'cm^-1')),
        HinderedRotor(inertia=(0.221535,'amu*angstrom^2'), symmetry=1, barrier=(13.5753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43223,0.0119348,-9.20111e-06,4.32805e-09,-8.54519e-13,24777.3,10.782], Tmin=(100,'K'), Tmax=(1197.14,'K')), NASAPolynomial(coeffs=[5.14053,0.00622677,-2.04895e-06,3.45073e-10,-2.27385e-14,24368.2,2.23311], Tmin=(1197.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CsJOOH)"""),
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
    label = '[O]O[CH]C=O(2043)',
    structure = SMILES('[O]O[CH]C=O'),
    E0 = (38.8345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.304072,'amu*angstrom^2'), symmetry=1, barrier=(6.99121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17369,'amu*angstrom^2'), symmetry=1, barrier=(56.7408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48922,0.0359335,-5.0258e-05,4.15082e-08,-1.40868e-11,4722.54,16.7287], Tmin=(100,'K'), Tmax=(782.765,'K')), NASAPolynomial(coeffs=[5.94543,0.0159796,-7.62826e-06,1.46012e-09,-1.01348e-13,4251.68,1.34957], Tmin=(782.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(OCJC=O)"""),
)

species(
    label = '[O]OC[C]=O(2044)',
    structure = SMILES('[O]OC[C]=O'),
    E0 = (55.0207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.654014,'amu*angstrom^2'), symmetry=1, barrier=(15.0371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340518,'amu*angstrom^2'), symmetry=1, barrier=(7.82917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63093,0.0324223,-3.97506e-05,2.80677e-08,-8.25459e-12,6664.69,17.1477], Tmin=(100,'K'), Tmax=(818.661,'K')), NASAPolynomial(coeffs=[6.29653,0.0145127,-6.93676e-06,1.34708e-09,-9.50121e-14,6064.49,0.19668], Tmin=(818.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.0207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C[CH]OO(2045)',
    structure = SMILES('O=C[CH]OO'),
    E0 = (-113.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42271,0.0370978,-3.65222e-05,1.95243e-08,-4.39492e-12,-13556.4,16.2149], Tmin=(100,'K'), Tmax=(1042.12,'K')), NASAPolynomial(coeffs=[7.46153,0.0177572,-8.68373e-06,1.71544e-09,-1.22643e-13,-14606.6,-8.30244], Tmin=(1042.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-113.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O)"""),
)

species(
    label = 'hoo_acetyl(2046)',
    structure = SMILES('O=[C]COO'),
    E0 = (-92.5524,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.855574,'amu*angstrom^2'), symmetry=1, barrier=(19.6713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855054,'amu*angstrom^2'), symmetry=1, barrier=(19.6594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855121,'amu*angstrom^2'), symmetry=1, barrier=(19.6609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93015,0.0435985,-5.03168e-05,2.8606e-08,-6.31892e-12,-11055.3,17.9195], Tmin=(100,'K'), Tmax=(1113.09,'K')), NASAPolynomial(coeffs=[11.5568,0.0090041,-3.69722e-06,6.84e-10,-4.76218e-14,-13198.3,-29.555], Tmin=(1113.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-92.5524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), label="""hoo_acetyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=COO[O](1387)',
    structure = SMILES('C=COO[O]'),
    E0 = (145.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,480.725,481.124,481.27,481.307,481.336],'cm^-1')),
        HinderedRotor(inertia=(0.000727217,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473985,'amu*angstrom^2'), symmetry=1, barrier=(7.79623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49568,0.0241603,5.63395e-06,-3.21506e-08,1.64367e-11,17602.2,19.5351], Tmin=(100,'K'), Tmax=(914.695,'K')), NASAPolynomial(coeffs=[12.482,0.00398125,2.01998e-07,-1.14149e-10,6.68075e-15,14792.6,-33.1245], Tmin=(914.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[O]OC=CO(2047)',
    structure = SMILES('[O]OC=CO'),
    E0 = (-98.9365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.971047,'amu*angstrom^2'), symmetry=1, barrier=(22.3263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968319,'amu*angstrom^2'), symmetry=1, barrier=(22.2635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647161,0.0499479,-5.65942e-05,2.8545e-08,-5.08564e-12,-11757.9,21.9324], Tmin=(100,'K'), Tmax=(1718.67,'K')), NASAPolynomial(coeffs=[15.3978,-0.00194453,4.02361e-06,-9.14063e-10,6.43635e-14,-14234.4,-49.673], Tmin=(1718.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.9365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(ROOJ)"""),
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
    E0 = (43.8946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (240.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (152.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (243.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (138.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (26.5684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-2.0894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (238.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (250.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (266.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (63.9354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-8.66325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (459.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (44.9303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction4',
    reactants = ['vinoxyperoxy(1385)'],
    products = ['HO2(9)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(2039)', 'CH3OO(44)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;C_pri] for rate rule [CO;C_pri/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', 'OCH2CHO(1527)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C[CH]O[O](2040)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['vinoxyperoxy(1385)'],
    products = ['[CH]1COOO1(2041)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(233.58,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 232.3 to 233.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['vinoxyperoxy(1385)'],
    products = ['[O]C1COO1(2042)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', 'vinoxy(1351)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(203262,'m^3/(mol*s)'), n=0.35323, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0148472730165, var=2.75207767881, Tref=1000.0, N=9, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R
    Total Standard Deviation in ln(k): 3.36303735057
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCO(1372)', '[CH2]O[O](46)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[O]O[CH]C=O(2043)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[O]OC[C]=O(2044)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['vinoxyperoxy(1385)'],
    products = ['O=C[CH]OO(2045)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_O;O_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['hoo_acetyl(2046)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;O_H_out] for rate rule [R4H_SSS;CO_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=COO[O](1387)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OC=CO(2047)'],
    products = ['vinoxyperoxy(1385)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

network(
    label = 'PDepNetwork #642',
    isomers = [
        'vinoxyperoxy(1385)',
    ],
    reactants = [
        ('HO2(9)', 'ketene(1375)'),
        ('O2(2)', 'vinoxy(1351)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #642',
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


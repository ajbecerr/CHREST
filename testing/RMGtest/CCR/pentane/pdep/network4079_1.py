species(
    label = 'C=[C][CH]O[C]=C[O](16660)',
    structure = SMILES('C=[C][CH]O[C]=C[O]'),
    E0 = (445.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1670,1700,300,440,408.891,408.892,408.895,408.899,408.9],'cm^-1')),
        HinderedRotor(inertia=(0.214908,'amu*angstrom^2'), symmetry=1, barrier=(25.4973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214903,'amu*angstrom^2'), symmetry=1, barrier=(25.4973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214901,'amu*angstrom^2'), symmetry=1, barrier=(25.4972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988273,0.0571091,-4.30181e-05,5.00785e-09,4.89842e-12,53667.6,28.403], Tmin=(100,'K'), Tmax=(956.338,'K')), NASAPolynomial(coeffs=[17.2899,0.01029,-3.09282e-06,5.35499e-10,-3.89945e-14,49572.6,-54.6231], Tmin=(956.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
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
    label = '[C]=C[O](6859)',
    structure = SMILES('[C]C=O'),
    E0 = (491.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53714,0.00730776,2.32451e-06,-8.52039e-09,3.75032e-12,59149.2,7.56143], Tmin=(100,'K'), Tmax=(1018.09,'K')), NASAPolynomial(coeffs=[6.47558,0.00262438,-8.8465e-07,2.00876e-10,-1.68063e-14,58195.3,-8.41394], Tmin=(1018.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJ3)"""),
)

species(
    label = '[CH2][C]=C[O](16852)',
    structure = SMILES('[CH2][C]=C[O]'),
    E0 = (363.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.0425533,'amu*angstrom^2'), symmetry=1, barrier=(25.7029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92133,0.0261555,-3.06573e-05,2.657e-08,-1.03537e-11,43797.7,13.4805], Tmin=(100,'K'), Tmax=(698.339,'K')), NASAPolynomial(coeffs=[4.09794,0.0176156,-8.44678e-06,1.67507e-09,-1.1987e-13,43677.3,8.54082], Tmin=(698.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH]=[C]O[CH][C]=C(16181)',
    structure = SMILES('[CH]=[C]O[CH][C]=C'),
    E0 = (759.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,180,597.029,833.203],'cm^-1')),
        HinderedRotor(inertia=(0.74415,'amu*angstrom^2'), symmetry=1, barrier=(17.1095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743464,'amu*angstrom^2'), symmetry=1, barrier=(17.0937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191273,'amu*angstrom^2'), symmetry=1, barrier=(95.2621,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3208.06,'J/mol'), sigma=(5.53191,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.09 K, Pc=43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69949,0.0499227,-5.16581e-05,2.80139e-08,-6.09932e-12,91451,24.9866], Tmin=(100,'K'), Tmax=(1109.68,'K')), NASAPolynomial(coeffs=[10.8034,0.0171054,-7.2963e-06,1.3616e-09,-9.46263e-14,89430.6,-19.8815], Tmin=(1109.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(759.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = 'C=[C][C]O[C]=C[O](18808)',
    structure = SMILES('[CH2][C]=[C]O[C]=C[O]'),
    E0 = (735.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,426.921,426.921,426.922,426.922],'cm^-1')),
        HinderedRotor(inertia=(0.924919,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131194,'amu*angstrom^2'), symmetry=1, barrier=(16.9684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131193,'amu*angstrom^2'), symmetry=1, barrier=(16.9684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08476,0.0601632,-6.9879e-05,4.01871e-08,-8.92995e-12,88543,31.5898], Tmin=(100,'K'), Tmax=(1111.78,'K')), NASAPolynomial(coeffs=[14.6599,0.0113217,-3.9821e-06,6.72429e-10,-4.44368e-14,85524.6,-35.3407], Tmin=(1111.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(735.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C1OC1=C[O](18809)',
    structure = SMILES('C=[C]C1OC1=C[O]'),
    E0 = (252.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11855,0.0414349,2.4012e-05,-8.15959e-08,3.99724e-11,30520,21.8457], Tmin=(100,'K'), Tmax=(923.326,'K')), NASAPolynomial(coeffs=[24.6577,-0.00362592,4.75438e-06,-9.31497e-10,5.58102e-14,23747,-102.977], Tmin=(923.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]OC1=CO1(18810)',
    structure = SMILES('C=C=CO[C]1[CH]O1'),
    E0 = (333.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12606,0.0770249,-8.49917e-05,4.30758e-08,-7.74352e-12,40311.1,29.7477], Tmin=(100,'K'), Tmax=(1709.17,'K')), NASAPolynomial(coeffs=[19.7256,0.00129168,5.11082e-06,-1.28886e-09,9.42914e-14,37117.3,-70.5175], Tmin=(1709.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Ethylene_oxide) + radical(Cs_P) + radical(CCsJO)"""),
)

species(
    label = 'C=C=COC=[C][O](16721)',
    structure = SMILES('C=C=CO[CH][C]=O'),
    E0 = (237.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16217,'amu*angstrom^2'), symmetry=1, barrier=(26.7205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15621,'amu*angstrom^2'), symmetry=1, barrier=(26.5835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15668,'amu*angstrom^2'), symmetry=1, barrier=(26.5943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.350094,0.0755191,-8.44133e-05,4.26662e-08,-7.91554e-12,28688.2,27.5296], Tmin=(100,'K'), Tmax=(1528.75,'K')), NASAPolynomial(coeffs=[23.5315,0.00025308,1.97622e-06,-4.75124e-10,3.36528e-14,22879.7,-92.9377], Tmin=(1528.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
)

species(
    label = 'C#C[CH]O[C]=C[O](18811)',
    structure = SMILES('[CH]=C=CO[C]=C[O]'),
    E0 = (437.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,540,610,2055,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.42208,'amu*angstrom^2'), symmetry=1, barrier=(32.6963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42859,'amu*angstrom^2'), symmetry=1, barrier=(32.8461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.237607,0.0654419,-7.3432e-05,3.83459e-08,-7.32445e-12,52743.4,30.4289], Tmin=(100,'K'), Tmax=(1516.12,'K')), NASAPolynomial(coeffs=[19.0071,0.00262992,1.86298e-06,-5.45422e-10,4.20748e-14,48579.7,-62.8962], Tmin=(1516.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C][CH]O[C]=C=O(18812)',
    structure = SMILES('C=[C][CH]O[C]=C=O'),
    E0 = (485.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,2120,512.5,787.5,384.414,384.415,384.416,384.417],'cm^-1')),
        HinderedRotor(inertia=(0.00114078,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16693,'amu*angstrom^2'), symmetry=1, barrier=(17.5054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646079,'amu*angstrom^2'), symmetry=1, barrier=(67.7508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46357,0.0584824,-7.4398e-05,5.06103e-08,-1.38241e-11,58452.3,26.1], Tmin=(100,'K'), Tmax=(892.142,'K')), NASAPolynomial(coeffs=[10.1451,0.0195572,-8.94998e-06,1.70235e-09,-1.18701e-13,56903.3,-14.7925], Tmin=(892.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'CH2CCHO(16849)',
    structure = SMILES('C=[C]C=O'),
    E0 = (174.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.41889,'amu*angstrom^2'), symmetry=1, barrier=(32.6232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.94221,0.0199513,-7.73714e-06,-1.79865e-09,1.45937e-12,21035.8,13.0042], Tmin=(100,'K'), Tmax=(1153.82,'K')), NASAPolynomial(coeffs=[6.85877,0.0121083,-4.99617e-06,9.2516e-10,-6.41112e-14,19750.3,-8.10553], Tmin=(1153.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""CH2CCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=[C]OC=C=C(16179)',
    structure = SMILES('C#CO[CH][C]=C'),
    E0 = (482.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100,279.596,279.598,279.634],'cm^-1')),
        HinderedRotor(inertia=(0.594372,'amu*angstrom^2'), symmetry=1, barrier=(32.9751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594396,'amu*angstrom^2'), symmetry=1, barrier=(32.9753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594461,'amu*angstrom^2'), symmetry=1, barrier=(32.9754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36759,0.0525792,-5.21763e-05,2.56965e-08,-4.95256e-12,58089.6,20.1027], Tmin=(100,'K'), Tmax=(1266.22,'K')), NASAPolynomial(coeffs=[13.9321,0.0128878,-5.15679e-06,9.40687e-10,-6.48368e-14,54907.7,-43.4797], Tmin=(1266.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH]O[C]=C[O](18813)',
    structure = SMILES('[CH]=[C][CH]O[C]=C[O]'),
    E0 = (692.343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,365.721,365.737,365.762,365.782],'cm^-1')),
        HinderedRotor(inertia=(0.279515,'amu*angstrom^2'), symmetry=1, barrier=(26.53,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279461,'amu*angstrom^2'), symmetry=1, barrier=(26.5303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279552,'amu*angstrom^2'), symmetry=1, barrier=(26.5303,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.603844,0.0642545,-7.21432e-05,3.8254e-08,-7.64793e-12,83400.9,30.1879], Tmin=(100,'K'), Tmax=(1332.48,'K')), NASAPolynomial(coeffs=[18.5824,0.00569444,-1.05413e-06,1.01559e-10,-4.75957e-15,79017.1,-60.18], Tmin=(1332.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(692.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH]O[C]=[C][O](18814)',
    structure = SMILES('C=[C][CH]O[C]=[C][O]'),
    E0 = (684.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,180,220.01,616.59,616.603,618.067],'cm^-1')),
        HinderedRotor(inertia=(0.0602668,'amu*angstrom^2'), symmetry=1, barrier=(16.2672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322548,'amu*angstrom^2'), symmetry=1, barrier=(87.1791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0601223,'amu*angstrom^2'), symmetry=1, barrier=(16.2684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31181,0.0574926,-6.49211e-05,3.70513e-08,-8.30614e-12,82483.7,30.3726], Tmin=(100,'K'), Tmax=(1091.57,'K')), NASAPolynomial(coeffs=[12.9932,0.0146863,-6.09721e-06,1.12465e-09,-7.78362e-14,79933.5,-27.0066], Tmin=(1091.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]O[C]=C[O](18815)',
    structure = SMILES('[CH]=C[CH]O[C]=C[O]'),
    E0 = (454.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,346.891,347.396,347.871,348.087],'cm^-1')),
        HinderedRotor(inertia=(0.314467,'amu*angstrom^2'), symmetry=1, barrier=(27.1625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316958,'amu*angstrom^2'), symmetry=1, barrier=(27.1681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316162,'amu*angstrom^2'), symmetry=1, barrier=(27.1665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950351,0.0559933,-3.47218e-05,-7.1505e-09,1.00946e-11,54783.8,28.4688], Tmin=(100,'K'), Tmax=(942.982,'K')), NASAPolynomial(coeffs=[18.7116,0.00797144,-1.78975e-06,2.89814e-10,-2.29759e-14,50219.5,-62.6165], Tmin=(942.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH]O[C]=[C]O(18816)',
    structure = SMILES('C=[C][CH]O[C]=[C]O'),
    E0 = (543.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,409.149,409.432,409.66,409.892],'cm^-1')),
        HinderedRotor(inertia=(0.140047,'amu*angstrom^2'), symmetry=1, barrier=(16.6655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140051,'amu*angstrom^2'), symmetry=1, barrier=(16.6645,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13993,'amu*angstrom^2'), symmetry=1, barrier=(16.664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00508,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.783638,0.065874,-7.34405e-05,3.77622e-08,-6.6619e-12,65491.6,30.8559], Tmin=(100,'K'), Tmax=(947.857,'K')), NASAPolynomial(coeffs=[16.4477,0.0106527,-3.27243e-06,5.22104e-10,-3.42632e-14,62033.3,-46.454], Tmin=(947.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][CH]OC=[C][O](16723)',
    structure = SMILES('C=[C][CH]OC=[C][O]'),
    E0 = (445.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,408.793,408.808,408.833,408.872,408.964],'cm^-1')),
        HinderedRotor(inertia=(0.214795,'amu*angstrom^2'), symmetry=1, barrier=(25.4961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214751,'amu*angstrom^2'), symmetry=1, barrier=(25.4965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21495,'amu*angstrom^2'), symmetry=1, barrier=(25.4992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988273,0.0571091,-4.30181e-05,5.00785e-09,4.89842e-12,53667.6,28.403], Tmin=(100,'K'), Tmax=(956.338,'K')), NASAPolynomial(coeffs=[17.2899,0.01029,-3.09282e-06,5.35499e-10,-3.89945e-14,49572.6,-54.6231], Tmin=(956.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]CO[C]=C[O](16661)',
    structure = SMILES('[CH]=[C]CO[C]=C[O]'),
    E0 = (581.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,207.085,207.105,207.116,207.135],'cm^-1')),
        HinderedRotor(inertia=(0.697925,'amu*angstrom^2'), symmetry=1, barrier=(21.2385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69785,'amu*angstrom^2'), symmetry=1, barrier=(21.2376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696898,'amu*angstrom^2'), symmetry=1, barrier=(21.238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.165229,0.0692637,-7.63312e-05,3.9347e-08,-7.54185e-12,70077.7,30.8289], Tmin=(100,'K'), Tmax=(1452.8,'K')), NASAPolynomial(coeffs=[20.3782,0.00435093,2.51764e-07,-1.83282e-10,1.55886e-14,65181.9,-70.8728], Tmin=(1452.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]CO[C]=[C][O](16662)',
    structure = SMILES('C=[C]CO[C]=[C][O]'),
    E0 = (574.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,424.171,424.195,424.217,424.235,424.25],'cm^-1')),
        HinderedRotor(inertia=(0.121122,'amu*angstrom^2'), symmetry=1, barrier=(15.4644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121102,'amu*angstrom^2'), symmetry=1, barrier=(15.4648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121113,'amu*angstrom^2'), symmetry=1, barrier=(15.4641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907683,0.0621482,-6.80727e-05,3.7006e-08,-7.77926e-12,69158.9,30.8857], Tmin=(100,'K'), Tmax=(1173.96,'K')), NASAPolynomial(coeffs=[15.5976,0.0120959,-4.12007e-06,6.88898e-10,-4.54218e-14,65709.8,-42.3408], Tmin=(1173.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=C[CH]O[C]=[C][O](18817)',
    structure = SMILES('C=C[CH]O[C]=[C][O]'),
    E0 = (447.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19815,0.0546398,-4.63161e-05,1.58925e-08,-8.37235e-13,53886.6,30.3053], Tmin=(100,'K'), Tmax=(1021.91,'K')), NASAPolynomial(coeffs=[14.5969,0.014485,-5.41674e-06,9.80793e-10,-6.86618e-14,50506.4,-37.7665], Tmin=(1021.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C][CH]OC=C[O](18818)',
    structure = SMILES('[CH]=[C][CH]OC=C[O]'),
    E0 = (452.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31906,'amu*angstrom^2'), symmetry=1, barrier=(30.3278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31997,'amu*angstrom^2'), symmetry=1, barrier=(30.3486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3236,'amu*angstrom^2'), symmetry=1, barrier=(30.4321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.714214,0.0587538,-3.23341e-05,-1.70569e-08,1.55226e-11,54565.9,26.6617], Tmin=(100,'K'), Tmax=(926.708,'K')), NASAPolynomial(coeffs=[21.6227,0.00341262,7.41094e-07,-2.03862e-10,1.06735e-14,49191.8,-80.7055], Tmin=(926.708,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]O[C]=CO(18819)',
    structure = SMILES('[CH]=[C][CH]O[C]=CO'),
    E0 = (550.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,201.949,201.949,201.949],'cm^-1')),
        HinderedRotor(inertia=(0.816452,'amu*angstrom^2'), symmetry=1, barrier=(23.6286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816451,'amu*angstrom^2'), symmetry=1, barrier=(23.6286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81645,'amu*angstrom^2'), symmetry=1, barrier=(23.6286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81645,'amu*angstrom^2'), symmetry=1, barrier=(23.6286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0574015,0.0742133,-8.62359e-05,4.63007e-08,-9.18405e-12,66414.6,31.1482], Tmin=(100,'K'), Tmax=(1414.01,'K')), NASAPolynomial(coeffs=[21.4907,0.00243952,1.37842e-06,-4.17545e-10,3.24399e-14,61402.3,-76.4496], Tmin=(1414.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=COC1=C[O](18820)',
    structure = SMILES('C=C1[CH]OC1=C[O]'),
    E0 = (86.5843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4976,0.0364381,2.69423e-05,-7.27464e-08,3.35647e-11,10520.9,19.1834], Tmin=(100,'K'), Tmax=(939.412,'K')), NASAPolynomial(coeffs=[19.5705,0.00656017,-5.18852e-07,8.62966e-11,-1.40013e-14,5048.05,-77.9343], Tmin=(939.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.5843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH2]C#CO[C]=C[O](18821)',
    structure = SMILES('[CH2]C#CO[C]=C[O]'),
    E0 = (465.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,2100,2250,500,550,371.23,371.354,371.379,371.457],'cm^-1')),
        HinderedRotor(inertia=(0.2835,'amu*angstrom^2'), symmetry=1, barrier=(27.7458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283457,'amu*angstrom^2'), symmetry=1, barrier=(27.7483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283562,'amu*angstrom^2'), symmetry=1, barrier=(27.7483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20906,0.0538174,-4.62324e-05,1.29089e-08,1.16404e-12,56127.7,27.3887], Tmin=(100,'K'), Tmax=(976.219,'K')), NASAPolynomial(coeffs=[15.9638,0.00959694,-3.23369e-06,5.81803e-10,-4.21047e-14,52473.3,-47.4019], Tmin=(976.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=COJ) + radical(Propargyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]O[C]=C[O](18822)',
    structure = SMILES('[CH2]C=[C]O[C]=C[O]'),
    E0 = (497.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,430.686,430.725,430.839,430.987],'cm^-1')),
        HinderedRotor(inertia=(0.136636,'amu*angstrom^2'), symmetry=1, barrier=(17.9879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136574,'amu*angstrom^2'), symmetry=1, barrier=(17.9936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234045,'amu*angstrom^2'), symmetry=1, barrier=(30.8171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13894,0.0553127,-4.42182e-05,9.80879e-09,2.49226e-12,59938.8,30.9219], Tmin=(100,'K'), Tmax=(959.355,'K')), NASAPolynomial(coeffs=[15.6638,0.0121335,-3.88288e-06,6.65342e-10,-4.65719e-14,56352,-42.7186], Tmin=(959.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C]OC=C[O](18823)',
    structure = SMILES('[CH2][C]=[C]OC=C[O]'),
    E0 = (495.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,274.589,274.589,274.589,274.589],'cm^-1')),
        HinderedRotor(inertia=(0.488483,'amu*angstrom^2'), symmetry=1, barrier=(26.1362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488483,'amu*angstrom^2'), symmetry=1, barrier=(26.1362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488483,'amu*angstrom^2'), symmetry=1, barrier=(26.1362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.225686,0.0658967,-6.85039e-05,3.36347e-08,-6.14445e-12,59750.7,31.5554], Tmin=(100,'K'), Tmax=(1541.31,'K')), NASAPolynomial(coeffs=[19.5297,0.00579036,-2.68102e-07,-9.25547e-11,9.44415e-15,54988.8,-66.0704], Tmin=(1541.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[C]=[C]O[C]=C[O](18824)',
    structure = SMILES('C[C]=[C]O[C]=C[O]'),
    E0 = (583.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,371.231,371.244,371.249,371.293],'cm^-1')),
        HinderedRotor(inertia=(0.134565,'amu*angstrom^2'), symmetry=1, barrier=(13.1566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134496,'amu*angstrom^2'), symmetry=1, barrier=(13.1568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134561,'amu*angstrom^2'), symmetry=1, barrier=(13.1566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2557,0.060595,-7.10276e-05,4.37474e-08,-1.07002e-11,70312.2,30.7], Tmin=(100,'K'), Tmax=(998.754,'K')), NASAPolynomial(coeffs=[11.8084,0.0183309,-7.55137e-06,1.37647e-09,-9.41442e-14,68204.3,-20.1975], Tmin=(998.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C]O[C]=CO(18825)',
    structure = SMILES('[CH2][C]=[C]O[C]=CO'),
    E0 = (593.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,305.943,305.943,305.943],'cm^-1')),
        HinderedRotor(inertia=(0.286631,'amu*angstrom^2'), symmetry=1, barrier=(19.0384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28663,'amu*angstrom^2'), symmetry=1, barrier=(19.0384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28663,'amu*angstrom^2'), symmetry=1, barrier=(19.0384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28663,'amu*angstrom^2'), symmetry=1, barrier=(19.0384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721131,0.066548,-7.11411e-05,3.10629e-08,-2.88554e-12,71543.9,31.4865], Tmin=(100,'K'), Tmax=(902.08,'K')), NASAPolynomial(coeffs=[17.6841,0.00802363,-1.58286e-06,1.70624e-10,-9.23183e-15,67804.3,-52.3666], Tmin=(902.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[C]=CO[C]=[C][O](18826)',
    structure = SMILES('C[C]=CO[C][C]=O'),
    E0 = (565.38,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.914316,'amu*angstrom^2'), symmetry=1, barrier=(21.0219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914165,'amu*angstrom^2'), symmetry=1, barrier=(21.0184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914307,'amu*angstrom^2'), symmetry=1, barrier=(21.0217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91428,'amu*angstrom^2'), symmetry=1, barrier=(21.0211,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156114,0.072123,-8.16727e-05,4.2821e-08,-8.43141e-12,68148.5,26.3845], Tmin=(100,'K'), Tmax=(1348.42,'K')), NASAPolynomial(coeffs=[21.426,0.00413623,-6.02271e-07,4.92327e-11,-2.69364e-15,62857,-80.9395], Tmin=(1348.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.38,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]C1=CO[C][CH]O1(18827)',
    structure = SMILES('C=C1[CH]O[C][CH]O1'),
    E0 = (471.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08961,0.0454105,6.1019e-06,-4.91161e-08,2.31206e-11,56788.7,18.7525], Tmin=(100,'K'), Tmax=(1005.51,'K')), NASAPolynomial(coeffs=[21.226,0.00893139,-4.55892e-06,1.10126e-09,-9.3222e-14,50533.9,-89.4708], Tmin=(1005.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(C=CCJ(O)C) + radical(CCsJOC(O)) + radical(CH2_triplet)"""),
)

species(
    label = '[C]1[CH]OC[C]=CO1(18828)',
    structure = SMILES('[C]1[CH]O[C]=COC1'),
    E0 = (572.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770994,0.00661154,0.000217164,-3.57744e-07,1.59543e-10,69053.5,23.127], Tmin=(100,'K'), Tmax=(891.158,'K')), NASAPolynomial(coeffs=[53.9941,-0.0600204,3.93672e-05,-7.82527e-09,5.27978e-13,52727.3,-265.889], Tmin=(891.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCsJOC(O)) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1=CO[C]C1[O](18829)',
    structure = SMILES('[CH2]C1=CO[C]C1[O]'),
    E0 = (502.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28764,0.0356591,4.50139e-05,-1.08268e-07,5.16272e-11,60552.6,18.2288], Tmin=(100,'K'), Tmax=(901.098,'K')), NASAPolynomial(coeffs=[25.1746,-0.00617353,7.77676e-06,-1.64998e-09,1.103e-13,53641.1,-108.987], Tmin=(901.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1[C]OC=[C]C1(18830)',
    structure = SMILES('[O]C1[C]OC=[C]C1'),
    E0 = (602.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49252,0.0373032,2.38654e-05,-7.33698e-08,3.57157e-11,72551,17.8731], Tmin=(100,'K'), Tmax=(909.737,'K')), NASAPolynomial(coeffs=[20.0124,0.0030928,2.41575e-06,-5.96874e-10,3.85214e-14,67227.4,-80.462], Tmin=(909.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CC(C)OJ) + radical(Cds_S) + radical(CH2_triplet)"""),
)

species(
    label = 'C=[C]C1O[C][CH]O1(18831)',
    structure = SMILES('C=[C]C1O[C][CH]O1'),
    E0 = (610.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.138456,0.0674185,-8.60683e-05,5.07427e-08,-1.05598e-11,73560.2,23.1629], Tmin=(100,'K'), Tmax=(1471.23,'K')), NASAPolynomial(coeffs=[16.661,-0.000287139,6.19074e-06,-1.58919e-09,1.21382e-13,71164.3,-54.5481], Tmin=(1471.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(1,3-Dioxolane) + radical(CCsJOCs) + radical(Cds_S) + radical(CH2_triplet)"""),
)

species(
    label = 'C=[C]C1O[C]C1[O](18832)',
    structure = SMILES('C=[C]C1O[C]C1[O]'),
    E0 = (727.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.21283,0.0593137,-5.66756e-05,2.68577e-08,-4.6665e-12,87672.5,27.2926], Tmin=(100,'K'), Tmax=(1725.36,'K')), NASAPolynomial(coeffs=[14.8084,0.00983809,-6.69891e-08,-2.68525e-10,2.51812e-14,84963.6,-44.3381], Tmin=(1725.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(Cds_S) + radical(CH2_triplet)"""),
)

species(
    label = 'C=C=CO[C]=C[O](16657)',
    structure = SMILES('C=C=CO[C]=C[O]'),
    E0 = (282.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25572,'amu*angstrom^2'), symmetry=1, barrier=(28.8714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2502,'amu*angstrom^2'), symmetry=1, barrier=(28.7446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982079,0.0555323,-3.43714e-05,-7.90855e-09,1.06823e-11,34132.4,27.1296], Tmin=(100,'K'), Tmax=(931.667,'K')), NASAPolynomial(coeffs=[18.6352,0.00733759,-1.20844e-06,1.54649e-10,-1.26313e-14,29645.4,-63.2146], Tmin=(931.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]O[C][CH][O](18833)',
    structure = SMILES('[CH2][CH][C]O[C]=C[O]'),
    E0 = (829.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.489045,0.0759897,-8.56588e-05,4.3325e-08,-7.97648e-12,99918.6,35.1596], Tmin=(100,'K'), Tmax=(1567.43,'K')), NASAPolynomial(coeffs=[23.6883,-0.0011884,3.01209e-06,-6.89267e-10,4.85377e-14,94240.7,-86.2831], Tmin=(1567.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(829.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJCO) + radical(RCCJ) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]O[C][CH][O](18834)',
    structure = SMILES('[CH][CH][CH]O[C]=C[O]'),
    E0 = (805.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.70079,0.0815138,-9.61908e-05,5.06483e-08,-9.68234e-12,97027.8,35.1765], Tmin=(100,'K'), Tmax=(1509.43,'K')), NASAPolynomial(coeffs=[24.8324,-0.00275624,4.05573e-06,-9.16201e-10,6.52636e-14,91211.6,-92.2526], Tmin=(1509.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(805.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CCJCO) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O][C][CH][O](10223)',
    structure = SMILES('[O][C][CH][O]'),
    E0 = (676.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,180,180,1781.81,1786.65],'cm^-1')),
        HinderedRotor(inertia=(0.371605,'amu*angstrom^2'), symmetry=1, barrier=(8.54393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61268,0.044967,-0.000114903,1.26746e-07,-4.85652e-11,81376.1,15.5028], Tmin=(100,'K'), Tmax=(883.837,'K')), NASAPolynomial(coeffs=[-0.0833967,0.0206941,-1.18058e-05,2.28838e-09,-1.54199e-13,83277.3,36.2363], Tmin=(883.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH) + radical(CH2_triplet)"""),
)

species(
    label = 'C=C=[C]O[C]C[O](18835)',
    structure = SMILES('[CH2]C#CO[C]C[O]'),
    E0 = (664.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2100,2250,500,550,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.911808,0.0701699,-8.70335e-05,5.43561e-08,-1.34327e-11,80079.8,24.0085], Tmin=(100,'K'), Tmax=(987.27,'K')), NASAPolynomial(coeffs=[13.5492,0.0189684,-9.24069e-06,1.82529e-09,-1.30633e-13,77584.5,-36.7977], Tmin=(987.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(664.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CCOJ) + radical(Propargyl) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C=CO[C]C[O](18836)',
    structure = SMILES('[CH]=C=CO[C]C[O]'),
    E0 = (636.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05393,'amu*angstrom^2'), symmetry=1, barrier=(24.2318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05795,'amu*angstrom^2'), symmetry=1, barrier=(24.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05921,'amu*angstrom^2'), symmetry=1, barrier=(24.3533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.309367,0.0776189,-0.000100364,6.26942e-08,-1.49633e-11,76679,25.712], Tmin=(100,'K'), Tmax=(1041.27,'K')), NASAPolynomial(coeffs=[17.8643,0.0101833,-3.22119e-06,4.99867e-10,-3.11635e-14,73023.1,-59.6904], Tmin=(1041.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(CH2_triplet) + radical(C=C=CJ)"""),
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
    E0 = (445.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (874.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (910.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1279.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (947.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (448.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (448.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (467.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (664.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (704.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (703.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (568.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (747.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (904.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (896.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (609.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (701.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (614.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (726.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (787.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (596.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (761.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (702.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (453.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (692.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (701.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (671.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (745.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (626.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (672.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (503.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (572.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (540.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (602.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (610.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (727.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (445.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (851.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (828.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (1014.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (814.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (787.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['OCHCO(3676)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O][C]=C[O](9592)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C[O](6859)', '[CH2][C]=C[O](16852)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]O[CH][C]=C(16181)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=[C][C]O[C]=C[O](18808)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=[C]C1OC1=C[O](18809)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=[C][CH]OC1=CO1(18810)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=C=COC=[C][O](16721)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C#C[CH]O[C]=C[O](18811)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C=[C][CH]O[C]=C=O(18812)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[C]=C[O](6859)', 'CH2CCHO(16849)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.86749e+08,'m^3/(mol*s)'), n=-0.730906, Ea=(37.2121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;Y_1centerbirad] + [Od_R;YJ] for rate rule [Od_R;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OCHCO(3676)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(29.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;Y_1centerbirad]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O(4)', '[CH]=[C]OC=C=C(16179)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.30864e+06,'m^3/(mol*s)'), n=-0.19959, Ea=(22.3126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_atom_triplet] + [Ct_Ct;YJ] for rate rule [Ct_Ct;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]=[C][CH]O[C]=C[O](18813)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C=[C][CH]O[C]=[C][O](18814)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C[CH]O[C]=C[O](18815)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][CH]O[C]=[C]O(18816)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]CO[C]=C[O](16661)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]CO[C]=[C][O](16662)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.65304e+08,'s^-1'), n=1.30038, Ea=(212.995,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_1;Cd_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=C[CH]O[C]=[C][O](18817)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5Hall;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C][CH]OC=C[O](18818)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C][CH]O[C]=CO(18819)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['[CH2]C1=COC1=C[O](18820)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]C#CO[C]=C[O](18821)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C=[C]O[C]=C[O](18822)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=[C]OC=C[O](18823)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_O;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[C]=[C]O[C]=C[O](18824)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=[C]O[C]=CO(18825)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[C]=CO[C]=[C][O](18826)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.48245e+06,'s^-1'), n=1.54936, Ea=(107.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out;Cs_H_out_2H] + [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['[CH2]C1=CO[C][CH]O1(18827)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(58.6284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['[C]1[CH]OC[C]=CO1(18828)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.97949e+09,'s^-1'), n=0.649948, Ea=(127.423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 122.5 to 127.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['[CH2]C1=CO[C]C1[O](18829)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.76201e+09,'s^-1'), n=0.626373, Ea=(94.9363,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['[O]C1[C]OC=[C]C1(18830)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.69238e+08,'s^-1'), n=0.607924, Ea=(157.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 154.7 to 157.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=[C]C1O[C][CH]O1(18831)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.25666e+09,'s^-1'), n=0.571459, Ea=(165.085,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 164.6 to 165.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=[C]C1O[C]C1[O](18832)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.32992e+07,'s^-1'), n=1.25825, Ea=(282.397,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;multiplebond_intra;radadd_intra_cs] for rate rule [R5;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 280.9 to 282.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C][CH]O[C]=C[O](16660)'],
    products = ['C=C=CO[C]=C[O](16657)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C=[C]O[C][CH][O](18833)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C[CH]O[C][CH][O](18834)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O][C][CH][O](10223)', 'C3H3(5450)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C=[C]O[C]C[O](18835)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.26683e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C=CO[C]C[O](18836)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4079',
    isomers = [
        'C=[C][CH]O[C]=C[O](16660)',
    ],
    reactants = [
        ('OCHCO(3676)', 'C3H3(5450)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4079',
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


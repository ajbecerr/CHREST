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
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C][C]OC=[C][O](18991)',
    structure = SMILES('[CH2][C]=[C]O[CH][C]=O'),
    E0 = (689.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1670,1700,300,440,1855,455,950,349.554,349.757,349.778],'cm^-1')),
        HinderedRotor(inertia=(0.233343,'amu*angstrom^2'), symmetry=1, barrier=(20.2838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233568,'amu*angstrom^2'), symmetry=1, barrier=(20.2848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233619,'amu*angstrom^2'), symmetry=1, barrier=(20.2818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233642,'amu*angstrom^2'), symmetry=1, barrier=(20.2892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499479,0.0715146,-9.04594e-05,5.34907e-08,-1.19906e-11,83065.9,29.2981], Tmin=(100,'K'), Tmax=(1110.96,'K')), NASAPolynomial(coeffs=[18.7348,0.00585804,-1.81064e-06,2.93936e-10,-1.96674e-14,79014.2,-60.5958], Tmin=(1110.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'C=[C]C1OC=C1[O](18862)',
    structure = SMILES('C=[C]C1OC=C1[O]'),
    E0 = (223.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23893,0.0383138,3.20505e-05,-8.92232e-08,4.25866e-11,27041.5,22.2756], Tmin=(100,'K'), Tmax=(921.751,'K')), NASAPolynomial(coeffs=[24.2028,-0.00311921,4.73196e-06,-9.40019e-10,5.66417e-14,20334.9,-100.057], Tmin=(921.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
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
    label = 'C#C[CH]OC=[C][O](18992)',
    structure = SMILES('[CH]=C=CO[CH][C]=O'),
    E0 = (391.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,540,610,2055,1855,455,950,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23723,'amu*angstrom^2'), symmetry=1, barrier=(28.4463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23848,'amu*angstrom^2'), symmetry=1, barrier=(28.4751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23808,'amu*angstrom^2'), symmetry=1, barrier=(28.4658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.309479,0.0763854,-9.27699e-05,5.0276e-08,-9.89248e-12,47264.5,27.997], Tmin=(100,'K'), Tmax=(1463.38,'K')), NASAPolynomial(coeffs=[23.1151,-0.00284509,4.02553e-06,-9.19486e-10,6.6347e-14,42036.5,-88.3702], Tmin=(1463.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCsJOC(O)) + radical(CsCJ=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C][O](6861)',
    structure = SMILES('[CH][C]=O'),
    E0 = (425.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1491.4,1492.35],'cm^-1')),
        HinderedRotor(inertia=(0.0736325,'amu*angstrom^2'), symmetry=1, barrier=(5.08901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62981,0.00805288,-4.35162e-06,1.03741e-09,-8.43183e-14,51134.2,10.4631], Tmin=(100,'K'), Tmax=(1941.7,'K')), NASAPolynomial(coeffs=[6.62042,0.00292938,-1.19496e-06,2.28726e-10,-1.56225e-14,49777.3,-6.4529], Tmin=(1941.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
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
    label = '[CH]=[C][CH]OC=[C][O](18993)',
    structure = SMILES('[CH][C]=CO[CH][C]=O'),
    E0 = (669.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0544591,0.0761913,-8.55948e-05,4.52086e-08,-8.98944e-12,80619.6,28.9884], Tmin=(100,'K'), Tmax=(1351,'K')), NASAPolynomial(coeffs=[21.3846,0.00642581,-1.15211e-06,9.35798e-11,-3.35649e-15,75400.7,-78.7684], Tmin=(1351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(Cds_S) + radical(CsCJ=O) + radical(AllylJ2_triplet)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3857.61,'J/mol'), sigma=(6.30207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.55 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988273,0.0571091,-4.30181e-05,5.00785e-09,4.89842e-12,53667.6,28.403], Tmin=(100,'K'), Tmax=(956.338,'K')), NASAPolynomial(coeffs=[17.2899,0.01029,-3.09282e-06,5.35499e-10,-3.89945e-14,49572.6,-54.6231], Tmin=(956.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]OC=[C][O](18994)',
    structure = SMILES('[CH]C=CO[CH][C]=O'),
    E0 = (431.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387852,0.0668053,-4.42593e-05,-5.25277e-09,1.09058e-11,51998.4,26.9252], Tmin=(100,'K'), Tmax=(939.64,'K')), NASAPolynomial(coeffs=[21.322,0.00903756,-2.08362e-06,3.28514e-10,-2.54662e-14,46680.4,-80.1305], Tmin=(939.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(AllylJ2_triplet) + radical(CsCJ=O)"""),
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
    label = '[CH]=[C]COC=[C][O](16724)',
    structure = SMILES('[CH]=[C]COC=[C][O]'),
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
    label = '[CH]=[C][CH]OC=[C]O(18995)',
    structure = SMILES('[CH]=[C][CH]OC=[C]O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0574015,0.0742133,-8.62359e-05,4.63007e-08,-9.18405e-12,66414.6,31.1482], Tmin=(100,'K'), Tmax=(1414.01,'K')), NASAPolynomial(coeffs=[21.4907,0.00243952,1.37842e-06,-4.17545e-10,3.24399e-14,61402.3,-76.4496], Tmin=(1414.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=COC1[C]=O(18996)',
    structure = SMILES('[CH2]C1=COC1[C]=O'),
    E0 = (153.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33378,0.0328406,5.3098e-05,-1.13792e-07,5.18483e-11,18629.7,22.1141], Tmin=(100,'K'), Tmax=(922.823,'K')), NASAPolynomial(coeffs=[25.4437,-0.00484397,5.7394e-06,-1.11445e-09,6.65551e-14,11334.7,-107.682], Tmin=(922.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C#CO[CH][C]=O(18997)',
    structure = SMILES('[CH2]C#CO[CH][C]=O'),
    E0 = (406.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2100,2250,500,550,1855,455,950,285.148,285.15,285.156],'cm^-1')),
        HinderedRotor(inertia=(0.459513,'amu*angstrom^2'), symmetry=1, barrier=(26.5132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459519,'amu*angstrom^2'), symmetry=1, barrier=(26.513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459482,'amu*angstrom^2'), symmetry=1, barrier=(26.5133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459468,'amu*angstrom^2'), symmetry=1, barrier=(26.513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242082,0.0679608,-7.66995e-05,3.9438e-08,-7.52588e-12,49048.6,26.4458], Tmin=(100,'K'), Tmax=(1443.1,'K')), NASAPolynomial(coeffs=[21.2654,0.00159458,6.96162e-07,-2.02665e-10,1.46588e-14,43823.6,-79.771], Tmin=(1443.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CCsJOCs) + radical(Propargyl) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]C=[C]O[CH][C]=O(18998)',
    structure = SMILES('[CH2]C=[C]O[CH][C]=O'),
    E0 = (451.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,1855,455,950,334.575,334.576,334.58],'cm^-1')),
        HinderedRotor(inertia=(0.268317,'amu*angstrom^2'), symmetry=1, barrier=(21.3148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268335,'amu*angstrom^2'), symmetry=1, barrier=(21.3148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268313,'amu*angstrom^2'), symmetry=1, barrier=(21.3149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268342,'amu*angstrom^2'), symmetry=1, barrier=(21.3149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.188161,0.0708922,-7.91969e-05,4.12232e-08,-8.04935e-12,54477.7,29.9467], Tmin=(100,'K'), Tmax=(1377.73,'K')), NASAPolynomial(coeffs=[20.934,0.00466767,-5.70815e-07,2.00112e-11,1.5994e-16,49330,-74.7238], Tmin=(1377.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(Allyl_P) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][C]=[C]OC[C]=O(18999)',
    structure = SMILES('[CH2][C]=[C]OC[C]=O'),
    E0 = (495.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,1855,455,950,298.42,298.422,298.425],'cm^-1')),
        HinderedRotor(inertia=(0.00189291,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387911,'amu*angstrom^2'), symmetry=1, barrier=(24.5146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387903,'amu*angstrom^2'), symmetry=1, barrier=(24.5145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387908,'amu*angstrom^2'), symmetry=1, barrier=(24.5146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787129,0.0630025,-6.63989e-05,3.40813e-08,-6.74825e-12,59733.2,28.7953], Tmin=(100,'K'), Tmax=(1245.15,'K')), NASAPolynomial(coeffs=[16.8921,0.0112658,-4.073e-06,7.11371e-10,-4.82605e-14,55722.6,-52.4331], Tmin=(1245.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = 'C[C]=[C]O[CH][C]=O(19000)',
    structure = SMILES('C[C]=[C]O[CH][C]=O'),
    E0 = (538.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,278.749,278.753,278.755],'cm^-1')),
        HinderedRotor(inertia=(0.284382,'amu*angstrom^2'), symmetry=1, barrier=(15.6813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284377,'amu*angstrom^2'), symmetry=1, barrier=(15.6812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284393,'amu*angstrom^2'), symmetry=1, barrier=(15.6814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284374,'amu*angstrom^2'), symmetry=1, barrier=(15.6813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686573,0.0717677,-9.1034e-05,5.63698e-08,-1.34938e-11,64834.4,28.3495], Tmin=(100,'K'), Tmax=(1031.16,'K')), NASAPolynomial(coeffs=[15.8627,0.0128966,-5.39438e-06,1.00095e-09,-6.9593e-14,61704.6,-45.3318], Tmin=(1031.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(Cds_S) + radical(C=CJO) + radical(CsCJ=O)"""),
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
    label = 'C=[C]C1OC1[C]=O(19001)',
    structure = SMILES('C=[C]C1OC1[C]=O'),
    E0 = (281.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72971,0.0366141,1.45106e-05,-5.98519e-08,3.12096e-11,33939.5,25.5262], Tmin=(100,'K'), Tmax=(874.285,'K')), NASAPolynomial(coeffs=[16.8232,0.00564275,2.30779e-06,-7.23463e-10,5.50373e-14,29844.8,-53.5873], Tmin=(874.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C][CH]OC[C]=O(19002)',
    structure = SMILES('[CH][C]=COC[C]=O'),
    E0 = (475.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712011,0.0621353,-4.26547e-05,2.08942e-09,6.0164e-12,57265.9,26.7614], Tmin=(100,'K'), Tmax=(966.759,'K')), NASAPolynomial(coeffs=[17.7107,0.0148808,-5.14314e-06,9.14281e-10,-6.50751e-14,52900.8,-60.2505], Tmin=(966.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CsCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=COC=C1[O](18973)',
    structure = SMILES('[CH2]C1=COC=C1[O]'),
    E0 = (-9.37808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71592,0.034562,2.03754e-05,-5.93962e-08,2.73725e-11,-1031.36,19.8155], Tmin=(100,'K'), Tmax=(945.934,'K')), NASAPolynomial(coeffs=[16.9374,0.0094847,-2.15994e-06,3.94215e-10,-3.39223e-14,-5668.81,-62.0643], Tmin=(945.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.37808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Furan) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
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
    label = '[CH2][C]=[C]OC=[C]O(19003)',
    structure = SMILES('[CH2][C]=[C]OC=[C]O'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721131,0.066548,-7.11411e-05,3.10629e-08,-2.88554e-12,71543.9,31.4865], Tmin=(100,'K'), Tmax=(902.08,'K')), NASAPolynomial(coeffs=[17.6841,0.00802363,-1.58286e-06,1.70624e-10,-9.23183e-15,67804.3,-52.3666], Tmin=(902.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH][C][CH2](17602)',
    structure = SMILES('[CH][C][CH2]'),
    E0 = (981.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,917.882,918.218,918.236],'cm^-1')),
        HinderedRotor(inertia=(0.00702614,'amu*angstrom^2'), symmetry=1, barrier=(4.20306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182792,'amu*angstrom^2'), symmetry=1, barrier=(4.20276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96474,0.0202917,-1.68078e-05,7.2211e-09,-1.21436e-12,118107,14.5039], Tmin=(100,'K'), Tmax=(1441.81,'K')), NASAPolynomial(coeffs=[8.09928,0.00604662,-1.9874e-06,3.68258e-10,-2.60955e-14,116627,-12.1458], Tmin=(1441.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(981.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]O[C]=C[O](19004)',
    structure = SMILES('[CH2][C][CH]O[C]=C[O]'),
    E0 = (810.516,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.540565,0.0838152,-0.000102928,5.6885e-08,-1.15381e-11,97659.8,33.4792], Tmin=(100,'K'), Tmax=(1383.78,'K')), NASAPolynomial(coeffs=[24.5993,-0.00139962,3.04239e-06,-7.19591e-10,5.25948e-14,91903.3,-91.6324], Tmin=(1383.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(810.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(RCCJ) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C1O[CH][C]1[O](18873)',
    structure = SMILES('C=[C]C1O[CH][C]1[O]'),
    E0 = (623.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.514949,0.06316,-6.89033e-05,3.76879e-08,-7.61388e-12,75155.6,27.007], Tmin=(100,'K'), Tmax=(1448.87,'K')), NASAPolynomial(coeffs=[15.1343,0.00962286,1.64588e-07,-3.69037e-10,3.58481e-14,72302.3,-44.1707], Tmin=(1448.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOCs) + radical(Cds_S)"""),
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
    label = '[O][C]1[CH]OC=[C]C1(18951)',
    structure = SMILES('[O]C1=CO[CH][C]C1'),
    E0 = (411.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23964,0.0404592,2.58818e-05,-8.58434e-08,4.34459e-11,49550.7,19.0419], Tmin=(100,'K'), Tmax=(890.615,'K')), NASAPolynomial(coeffs=[23.4867,-0.0038628,6.89441e-06,-1.53965e-09,1.06562e-13,43383.1,-98.0887], Tmin=(890.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(C=C(C)OJ) + radical(CCsJOC(O)) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]1[CH]OC=[C]CO1(18903)',
    structure = SMILES('[C]1[CH]OC=[C]OC1'),
    E0 = (572.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770994,0.00661154,0.000217164,-3.57744e-07,1.59543e-10,69053.5,23.127], Tmin=(100,'K'), Tmax=(891.158,'K')), NASAPolynomial(coeffs=[53.9941,-0.0600204,3.93672e-05,-7.82527e-09,5.27978e-13,52727.3,-265.889], Tmin=(891.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCsJOC(O)) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]C1OC1[C]=O(19005)',
    structure = SMILES('[CH2][C]C1OC1[C]=O'),
    E0 = (569.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234055,0.0604083,-6.03195e-05,2.99225e-08,-5.41713e-12,68659.3,30.733], Tmin=(100,'K'), Tmax=(1663.24,'K')), NASAPolynomial(coeffs=[14.682,0.0092465,6.25219e-07,-4.39438e-10,3.84505e-14,66123.8,-39.4949], Tmin=(1663.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C1C[C][CH]O1(19006)',
    structure = SMILES('O=[C]C1C[C][CH]O1'),
    E0 = (473.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82908,0.0302087,4.06636e-05,-9.18882e-08,4.39498e-11,57036.1,21.8711], Tmin=(100,'K'), Tmax=(877.651,'K')), NASAPolynomial(coeffs=[18.5828,0.00259419,4.55355e-06,-1.1791e-09,8.58891e-14,52218.1,-67.4642], Tmin=(877.651,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CCsJOCs) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C][CH][O](17607)',
    structure = SMILES('[CH2][C][CH][O]'),
    E0 = (786.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1272.15,2748.05,2748.59],'cm^-1')),
        HinderedRotor(inertia=(1.68583,'amu*angstrom^2'), symmetry=1, barrier=(38.7606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68553,'amu*angstrom^2'), symmetry=1, barrier=(38.7536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71728,0.0335269,-5.7437e-05,5.8545e-08,-2.28758e-11,94644.4,17.6454], Tmin=(100,'K'), Tmax=(814.316,'K')), NASAPolynomial(coeffs=[3.31265,0.0190527,-9.49991e-06,1.88219e-09,-1.32741e-13,94930.4,17.2466], Tmin=(814.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(786.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(RCCJ) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C][O](10218)',
    structure = SMILES('[CH][C][O]'),
    E0 = (885.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([418.405,418.464,1386.35,1386.4,1890.01],'cm^-1')),
        HinderedRotor(inertia=(0.0607177,'amu*angstrom^2'), symmetry=1, barrier=(7.54278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.35351,0.0168371,-3.18729e-05,3.43969e-08,-1.43789e-11,106493,11.3039], Tmin=(100,'K'), Tmax=(740.605,'K')), NASAPolynomial(coeffs=[3.99585,0.00842558,-4.82643e-06,1.03992e-09,-7.72054e-14,106534,9.313], Tmin=(740.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(885.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCJ2_triplet) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C=[C]O[CH][C][O](19007)',
    structure = SMILES('[CH2][CH][C]OC=[C][O]'),
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
    label = '[CH]=C[CH]O[CH][C][O](19008)',
    structure = SMILES('[CH][CH][CH]O[CH][C]=O'),
    E0 = (799.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.45718,0.0785991,-0.000108152,7.25929e-08,-1.88071e-11,96339.8,31.2306], Tmin=(100,'K'), Tmax=(954.742,'K')), NASAPolynomial(coeffs=[16.024,0.0133793,-5.68369e-06,1.04185e-09,-7.12243e-14,93367.3,-43.1493], Tmin=(954.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(799.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCJCO) + radical(CCsJOCs) + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C1[CH]OC1[C][O](19009)',
    structure = SMILES('[CH2]C1=COC1[C][O]'),
    E0 = (601.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.906812,0.0491404,9.78054e-07,-5.36251e-08,2.8447e-11,72510.1,22.2131], Tmin=(100,'K'), Tmax=(941.185,'K')), NASAPolynomial(coeffs=[23.3459,0.0013275,1.39302e-06,-2.37655e-10,7.08112e-15,66180.1,-95.8714], Tmin=(941.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CCOJ) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C1=CO[CH][C]O1(18921)',
    structure = SMILES('C=C1[CH]O[CH][C]O1'),
    E0 = (457.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1604,0.0429937,1.38452e-05,-5.79276e-08,2.65366e-11,55166.7,18.4694], Tmin=(100,'K'), Tmax=(994.145,'K')), NASAPolynomial(coeffs=[21.3031,0.00836496,-3.941e-06,9.64737e-10,-8.35593e-14,48868,-90.1256], Tmin=(994.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(C=CCJ(O)C) + radical(CCsJOCs) + radical(CH2_triplet)"""),
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
    E0 = (445.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (445.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (818.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (901.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (453.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (712.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (618.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (630.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (788.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (896.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (880.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (614.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (640.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (727.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (726.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (693.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (595.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (604.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (702.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (453.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (633.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (570.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (655.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (649.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (699.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (708.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (448.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (508.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (452.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (589.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (575.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (598.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (626.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (906.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (832.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (623.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (610.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (525.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (572.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (569.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (516.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (951.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (1059.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (851.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (822.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (601.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (508.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (1014.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['OCHCO(3676)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=C=COC=[C][O](16721)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['HCCO(2227)', 'CH2CCHO(16849)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O][C]=C[O](9592)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=[C][C]OC=[C][O](18991)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=[C]C1OC=C1[O](18862)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C=[C][CH]O[C]=C=O(18812)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C#C[CH]OC=[C][O](18992)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C][O](6861)', 'CH2CCHO(16849)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C][O](6861)', '[CH2][C]=C[O](16852)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'C=[C][CH]O[C]=[C][O](18814)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]=[C][CH]OC=[C][O](18993)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=[C][CH]O[C]=C[O](16660)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[CH]=C[CH]OC=[C][O](18994)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]CO[C]=[C][O](16662)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.67748e+09,'s^-1'), n=1.16185, Ea=(153.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]COC=[C][O](16724)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C][CH]O[C]=[C]O(18816)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=C[CH]O[C]=[C][O](18817)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;XH_out] for rate rule [R4HJ_1;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][CH]OC=C[O](18818)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C][CH]OC=[C]O(18995)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[CH2]C1=COC1[C]=O(18996)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C#CO[CH][C]=O(18997)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['OCHCO(3676)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C=[C]O[CH][C]=O(18998)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=[C]OC[C]=O(18999)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.67748e+09,'s^-1'), n=1.16185, Ea=(153.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[C]=[C]O[CH][C]=O(19000)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=[C]OC=C[O](18823)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.82652e+08,'s^-1'), n=1.30038, Ea=(212.995,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=[C]C1OC1[C]=O(19001)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C][CH]OC[C]=O(19002)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[CH2]C1=COC=C1[O](18973)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_DSSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O][C]=C[O](9592)', 'C3H3(5450)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(260000,'m^3/(mol*s)'), n=0, Ea=(46.5868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Ct_Ct;OJ_sec] for rate rule [Ct_Ct;O_rad/OneDe]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HCCO(2227)', '[CH2][C]=C[O](16852)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(260000,'m^3/(mol*s)'), n=0, Ea=(46.5868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Ct_Ct;OJ_sec] for rate rule [Ct_Ct;O_rad/OneDe]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[C]=CO[C]=[C][O](18826)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C]=[C]OC=[C]O(19003)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_3;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['OCHCO(3676)', '[CH][C][CH2](17602)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C][CH]O[C]=C[O](19004)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=[C]C1O[CH][C]1[O](18873)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(178.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 178.3 to 178.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=[C]C1O[C][CH]O1(18831)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.25666e+09,'s^-1'), n=0.571459, Ea=(165.085,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;multiplebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 164.6 to 165.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[O][C]1[CH]OC=[C]C1(18951)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.26172e+08,'s^-1'), n=0.58655, Ea=(80.0039,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[C]1[CH]OC=[C]CO1(18903)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.42391e+09,'s^-1'), n=0.475668, Ea=(127.423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7;multiplebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 122.5 to 127.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[CH2][C]C1OC1[C]=O(19005)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.86885e+09,'s^-1'), n=0.611527, Ea=(124.33,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 123.1 to 124.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['O=[C]C1C[C][CH]O1(19006)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.55592e+08,'s^-1'), n=0.712397, Ea=(71.4202,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['HCCO(2227)', '[CH2][C][CH][O](17607)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH][C][O](10218)', 'CH2CCHO(16849)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=[C]O[CH][C][O](19007)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C[CH]O[CH][C][O](19008)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['C=C1[CH]OC1[C][O](19009)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(156.566,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 155.9 to 156.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C][CH]OC=[C][O](16723)'],
    products = ['[CH2]C1=CO[CH][C]O1(18921)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O][C][CH][O](10223)', 'C3H3(5450)'],
    products = ['C=[C][CH]OC=[C][O](16723)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

network(
    label = 'PDepNetwork #4064',
    isomers = [
        'C=[C][CH]OC=[C][O](16723)',
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
    label = 'PDepNetwork #4064',
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


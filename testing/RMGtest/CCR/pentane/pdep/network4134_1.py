species(
    label = '[CH2][C]=CC=[C][CH2](16912)',
    structure = SMILES('[CH2][C]=CC=[C][CH2]'),
    E0 = (733.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,383.297],'cm^-1')),
        HinderedRotor(inertia=(2.93244,'amu*angstrom^2'), symmetry=1, barrier=(67.4226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191132,'amu*angstrom^2'), symmetry=1, barrier=(67.4298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188775,'amu*angstrom^2'), symmetry=1, barrier=(67.4008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84928,0.0398838,-8.2793e-06,-1.72322e-08,9.68702e-12,88323.6,22.6963], Tmin=(100,'K'), Tmax=(960.212,'K')), NASAPolynomial(coeffs=[10.3621,0.0219425,-7.62276e-06,1.3151e-09,-8.95683e-14,85881.1,-22.2335], Tmin=(960.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
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
    label = '[CH][C]=CC=[C][CH2](17824)',
    structure = SMILES('[CH][C]=CC=[C][CH2]'),
    E0 = (986.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,278.621,278.621,278.622,278.626],'cm^-1')),
        HinderedRotor(inertia=(0.924873,'amu*angstrom^2'), symmetry=1, barrier=(50.9482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92488,'amu*angstrom^2'), symmetry=1, barrier=(50.9482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.924852,'amu*angstrom^2'), symmetry=1, barrier=(50.9481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41779,0.0503908,-3.64914e-05,1.42455e-08,-2.30026e-12,118722,24.4952], Tmin=(100,'K'), Tmax=(1446.83,'K')), NASAPolynomial(coeffs=[11.0726,0.0236987,-8.81848e-06,1.49455e-09,-9.70116e-14,115928,-25.6497], Tmin=(1446.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(986.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC=C1[CH2](17825)',
    structure = SMILES('[CH2]C1=C[CH]C1=C'),
    E0 = (435.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14032,0.0143017,9.81667e-05,-1.50541e-07,6.26012e-11,52453.1,15.047], Tmin=(100,'K'), Tmax=(921.859,'K')), NASAPolynomial(coeffs=[20.1991,0.00455366,2.38964e-06,-5.4289e-10,2.86312e-14,46208.2,-86.4189], Tmin=(921.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC#C[CH2](17826)',
    structure = SMILES('[CH2][C]=CC#C[CH2]'),
    E0 = (708.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2100,2250,500,550,322.101],'cm^-1')),
        HinderedRotor(inertia=(0.00162484,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664835,'amu*angstrom^2'), symmetry=1, barrier=(48.946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664825,'amu*angstrom^2'), symmetry=1, barrier=(48.9461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87382,0.0430323,-3.17318e-05,1.24373e-08,-2.01431e-12,85252.2,20.9301], Tmin=(100,'K'), Tmax=(1434.79,'K')), NASAPolynomial(coeffs=[10.0068,0.0203588,-8.02804e-06,1.42354e-09,-9.52669e-14,82918.4,-21.243], Tmin=(1434.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CTCC=CCJ) + radical(Propargyl) + radical(Cds_S)"""),
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
    label = '[CH2][C]=[C]C=[C][CH2](17827)',
    structure = SMILES('[CH2][C]=[C]C=[C][CH2]'),
    E0 = (932.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1685,1700,300,370,440,436.956],'cm^-1')),
        HinderedRotor(inertia=(0.427222,'amu*angstrom^2'), symmetry=1, barrier=(57.8844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42723,'amu*angstrom^2'), symmetry=1, barrier=(57.8845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427245,'amu*angstrom^2'), symmetry=1, barrier=(57.8844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80503,0.0460057,-3.57856e-05,1.3417e-08,-1.05105e-12,112254,23.1467], Tmin=(100,'K'), Tmax=(899.621,'K')), NASAPolynomial(coeffs=[8.74056,0.0221725,-7.72557e-06,1.27763e-09,-8.26504e-14,110723,-11.1549], Tmin=(899.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(932.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]=C[CH2](17828)',
    structure = SMILES('[CH2][C]=C[C]=C[CH2]'),
    E0 = (694.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81618,0.0416734,-1.19731e-05,-1.45409e-08,9.33871e-12,83651.9,22.6322], Tmin=(100,'K'), Tmax=(924.203,'K')), NASAPolynomial(coeffs=[9.99496,0.02256,-7.38222e-06,1.21296e-09,-8.00138e-14,81444.6,-19.9437], Tmin=(924.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C=C[CH2](17829)',
    structure = SMILES('[CH2][C]=[C]C=C[CH2]'),
    E0 = (694.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81618,0.0416734,-1.19731e-05,-1.45409e-08,9.33871e-12,83651.9,22.6322], Tmin=(100,'K'), Tmax=(924.203,'K')), NASAPolynomial(coeffs=[9.99496,0.02256,-7.38222e-06,1.21296e-09,-8.00138e-14,81444.6,-19.9437], Tmin=(924.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]=[C]C(17830)',
    structure = SMILES('[CH2][C]=C[C]=[C]C'),
    E0 = (814.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180],'cm^-1')),
        HinderedRotor(inertia=(2.84496,'amu*angstrom^2'), symmetry=1, barrier=(65.4113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84143,'amu*angstrom^2'), symmetry=1, barrier=(65.3302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636631,'amu*angstrom^2'), symmetry=1, barrier=(14.6374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59009,0.0541157,-5.61868e-05,3.56219e-08,-9.48191e-12,98060.1,22.3624], Tmin=(100,'K'), Tmax=(902.277,'K')), NASAPolynomial(coeffs=[7.79607,0.0266062,-1.04583e-05,1.83806e-09,-1.22224e-13,96940.1,-6.94034], Tmin=(902.277,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=[C]C=[C]C(17831)',
    structure = SMILES('[CH2][C]=[C]C=[C]C'),
    E0 = (814.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180],'cm^-1')),
        HinderedRotor(inertia=(2.84496,'amu*angstrom^2'), symmetry=1, barrier=(65.4113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84143,'amu*angstrom^2'), symmetry=1, barrier=(65.3302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636631,'amu*angstrom^2'), symmetry=1, barrier=(14.6374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59009,0.0541157,-5.61868e-05,3.56219e-08,-9.48191e-12,98060.1,22.3624], Tmin=(100,'K'), Tmax=(902.277,'K')), NASAPolynomial(coeffs=[7.79607,0.0266062,-1.04583e-05,1.83806e-09,-1.22224e-13,96940.1,-6.94034], Tmin=(902.277,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(814.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC1[C]=C(17832)',
    structure = SMILES('C=[C]C1[CH]C1=C'),
    E0 = (637.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85862,0.0390446,-8.50405e-06,-1.45968e-08,7.71885e-12,76769.7,16.8694], Tmin=(100,'K'), Tmax=(1031.23,'K')), NASAPolynomial(coeffs=[10.669,0.0220019,-8.63336e-06,1.59647e-09,-1.12343e-13,74041.7,-30.3229], Tmin=(1031.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C=C[C]=C(17833)',
    structure = SMILES('[CH2]C=[C]C=C=C'),
    E0 = (515.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60107,0.0457135,-1.98528e-05,-8.88689e-09,7.69602e-12,62096.2,20.4258], Tmin=(100,'K'), Tmax=(942.99,'K')), NASAPolynomial(coeffs=[11.6413,0.0202885,-6.7122e-06,1.12524e-09,-7.56139e-14,59439.5,-31.4698], Tmin=(942.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[CH]C=[C][CH2](17834)',
    structure = SMILES('[CH]=C=CC=[C][CH2]'),
    E0 = (708.849,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.83763,'amu*angstrom^2'), symmetry=1, barrier=(42.2506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83615,'amu*angstrom^2'), symmetry=1, barrier=(42.2166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11113,0.051108,-4.51431e-05,2.08071e-08,-3.69076e-12,85369.4,23.6956], Tmin=(100,'K'), Tmax=(1566.09,'K')), NASAPolynomial(coeffs=[13.2742,0.0137116,-3.26162e-06,3.976e-10,-2.06987e-14,82336.1,-37.9617], Tmin=(1566.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2][C]=[C]C[C]=C(17471)',
    structure = SMILES('[CH2][C]=[C]C[C]=C'),
    E0 = (917.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,293.166,294.973],'cm^-1')),
        HinderedRotor(inertia=(0.00196266,'amu*angstrom^2'), symmetry=1, barrier=(0.119691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157357,'amu*angstrom^2'), symmetry=1, barrier=(9.60911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367991,'amu*angstrom^2'), symmetry=1, barrier=(22.4917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8074,0.0510772,-4.84015e-05,2.70682e-08,-6.53348e-12,110456,24.2111], Tmin=(100,'K'), Tmax=(969.403,'K')), NASAPolynomial(coeffs=[7.43063,0.027875,-1.25007e-05,2.37963e-09,-1.66713e-13,109366,-2.74314], Tmin=(969.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(917.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
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
    label = '[CH]=[C]CC=[C][CH2](16907)',
    structure = SMILES('[CH]=[C]CC=[C][CH2]'),
    E0 = (927.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0544244,'amu*angstrom^2'), symmetry=1, barrier=(10.4638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135165,'amu*angstrom^2'), symmetry=1, barrier=(25.0247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.348043,'amu*angstrom^2'), symmetry=1, barrier=(10.5051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68217,0.0511066,-4.45835e-05,2.12564e-08,-4.22835e-12,111576,24.5816], Tmin=(100,'K'), Tmax=(1179.05,'K')), NASAPolynomial(coeffs=[9.59155,0.0242736,-1.04464e-05,1.95443e-09,-1.35673e-13,109711,-14.8793], Tmin=(1179.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(927.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C=C[CH2](17836)',
    structure = SMILES('[CH][C]=CC=C[CH2]'),
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
    label = '[CH]=[C][CH]C=[C]C(17837)',
    structure = SMILES('[CH]=C=C[CH][C]C'),
    E0 = (844.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,337.04,338.827],'cm^-1')),
        HinderedRotor(inertia=(0.217906,'amu*angstrom^2'), symmetry=1, barrier=(17.8537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978072,'amu*angstrom^2'), symmetry=1, barrier=(80.9381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975377,'amu*angstrom^2'), symmetry=1, barrier=(80.9065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43001,0.0506075,-3.64345e-05,1.01877e-08,1.9734e-13,101617,21.8315], Tmin=(100,'K'), Tmax=(1019.56,'K')), NASAPolynomial(coeffs=[11.7668,0.020495,-7.49403e-06,1.30903e-09,-8.85985e-14,98966.7,-30.8992], Tmin=(1019.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(844.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH][C]=CC[C]=C(17468)',
    structure = SMILES('[CH][C]=CC[C]=C'),
    E0 = (899.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,382.104,382.105,382.105,382.105,382.106],'cm^-1')),
        HinderedRotor(inertia=(0.499374,'amu*angstrom^2'), symmetry=1, barrier=(51.7393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499379,'amu*angstrom^2'), symmetry=1, barrier=(51.7393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499376,'amu*angstrom^2'), symmetry=1, barrier=(51.7393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88044,0.0486308,-3.16095e-05,1.10825e-08,-1.71753e-12,108210,24.2874], Tmin=(100,'K'), Tmax=(1376.5,'K')), NASAPolynomial(coeffs=[7.4056,0.0325752,-1.41134e-05,2.60888e-09,-1.78548e-13,106689,-4.13387], Tmin=(1376.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(899.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=C[C][CH2](17838)',
    structure = SMILES('[CH]C=[C][CH2]'),
    E0 = (731.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,296.375,297.489,297.589],'cm^-1')),
        HinderedRotor(inertia=(0.828676,'amu*angstrom^2'), symmetry=1, barrier=(51.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815826,'amu*angstrom^2'), symmetry=1, barrier=(51.0347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,88067.9,15.929], Tmin=(100,'K'), Tmax=(2065.82,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,84627.8,-30.2813], Tmin=(2065.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(731.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]C1[CH]C1=C(17839)',
    structure = SMILES('[CH2][C]C1[CH]C1=C'),
    E0 = (925.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84146,0.0340958,1.84408e-05,-5.09924e-08,2.26163e-11,111460,20.0596], Tmin=(100,'K'), Tmax=(962.882,'K')), NASAPolynomial(coeffs=[13.9136,0.0167422,-5.61655e-06,1.03792e-09,-7.70593e-14,107615,-45.6197], Tmin=(962.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(925.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=[C]C1[CH][C]C1(17840)',
    structure = SMILES('C=[C]C1[CH][C]C1'),
    E0 = (960.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26864,0.0262024,2.87497e-05,-5.24121e-08,2.07103e-11,115536,22.8311], Tmin=(100,'K'), Tmax=(997.407,'K')), NASAPolynomial(coeffs=[10.4631,0.0222248,-8.70893e-06,1.66094e-09,-1.20888e-13,112465,-23.8835], Tmin=(997.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(960.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C1[C]C1(17841)',
    structure = SMILES('[CH2][C]=CC1[C]C1'),
    E0 = (915.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87564,0.0322298,2.50516e-05,-5.92592e-08,2.59435e-11,110241,21.7213], Tmin=(100,'K'), Tmax=(956.03,'K')), NASAPolynomial(coeffs=[14.5083,0.0153168,-4.80422e-06,8.84068e-10,-6.6974e-14,106183,-47.2465], Tmin=(956.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(915.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(Cds_S)"""),
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
    label = '[CH2][C][CH][C]=C[CH2](17842)',
    structure = SMILES('[CH2][C][CH][C]=C[CH2]'),
    E0 = (1107.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,308.572,1138.78,1139.09],'cm^-1')),
        HinderedRotor(inertia=(0.085833,'amu*angstrom^2'), symmetry=1, barrier=(5.58529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628097,'amu*angstrom^2'), symmetry=1, barrier=(42.4189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38114,'amu*angstrom^2'), symmetry=1, barrier=(93.9218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0272068,'amu*angstrom^2'), symmetry=1, barrier=(25.0321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4882,0.0514989,-4.12595e-05,1.71647e-08,-2.91243e-12,133305,25.1789], Tmin=(100,'K'), Tmax=(1385.39,'K')), NASAPolynomial(coeffs=[11.8013,0.021722,-9.01898e-06,1.65007e-09,-1.1272e-13,130448,-27.9376], Tmin=(1385.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1107.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P) + radical(RCCJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH][CH][C][CH2](17843)',
    structure = SMILES('[CH]C=C[CH][C][CH2]'),
    E0 = (1088.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23829,0.0523455,-3.38912e-05,1.09798e-08,-1.44397e-12,131075,26.4523], Tmin=(100,'K'), Tmax=(1750.43,'K')), NASAPolynomial(coeffs=[13.825,0.0235826,-9.24315e-06,1.59233e-09,-1.03219e-13,126668,-41.3182], Tmin=(1750.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1088.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1[CH][CH][C]C1(17844)',
    structure = SMILES('[CH2]C1=C[CH][C]C1'),
    E0 = (709.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2823,0.0182751,6.84229e-05,-1.03933e-07,4.13966e-11,85468.2,16.6665], Tmin=(100,'K'), Tmax=(956.338,'K')), NASAPolynomial(coeffs=[14.6149,0.0158622,-4.91417e-06,9.52444e-10,-7.59822e-14,80860.9,-54.0361], Tmin=(956.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-allyl) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]1[CH]C=[C]CC1(17845)',
    structure = SMILES('[C]1[CH]C=[C]CC1'),
    E0 = (803.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50773,0.0160034,6.50288e-05,-9.28772e-08,3.55523e-11,96694.3,15.1454], Tmin=(100,'K'), Tmax=(974.049,'K')), NASAPolynomial(coeffs=[11.6043,0.020947,-7.72327e-06,1.49935e-09,-1.13188e-13,92915.6,-38.8013], Tmin=(974.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(cyclohexene-allyl) + radical(CCJ2_triplet) + radical(Cds_S)"""),
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
    E0 = (733.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1198.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (741.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (935.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (975.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1228.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1144.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (883.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (927.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1007.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (983.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (736.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (756.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (935.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1100.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (903.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1072.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (781.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (976.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1108.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1387.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (925.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (960.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (915.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1375.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1130.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1111.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (850.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (846.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['C3H3(5450)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[CH][C]=CC=[C][CH2](17824)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2]C1=CC=C1[CH2](17825)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2][C]=CC#C[CH2](17826)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C3H3(5450)', '[CH]=[C][CH2](16918)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C][CH2](16918)', '[CH]=[C][CH2](16918)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.625e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C]=[C]C=[C][CH2](17827)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2][C]=C[C]=C[CH2](17828)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.68855e+09,'s^-1'), n=1.26608, Ea=(149.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2][C]=[C]C=C[CH2](17829)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=C[C]=[C]C(17830)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=[C]C=[C]C(17831)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(115.297,'s^-1'), n=2.99825, Ea=(169.04,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2]C1=CC1[C]=C(17832)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.36e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2]C=C=C[C]=C(17833)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'C#C[CH]C=[C][CH2](17834)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C]=[C]C[C]=C(17471)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.28974e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C[CH]C=[C][CH2](17835)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]CC=[C][CH2](16907)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C][CH]C=C[CH2](17836)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][CH]C=[C]C(17837)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH][C]=CC[C]=C(17468)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H2CC(T)(1341)', '[CH]=C[C][CH2](17838)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2][C]C1[CH]C1=C(17839)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.121e+13,'s^-1'), n=0.275, Ea=(192.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 191.7 to 192.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['C=[C]C1[CH][C]C1(17840)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.28131e+09,'s^-1'), n=0.807397, Ea=(226.353,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 224.4 to 226.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['C=[C][CH]C1[C]C1(17841)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.01734e+09,'s^-1'), n=0.788889, Ea=(182.188,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 181.3 to 182.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H3(5450)', '[CH][C][CH2](17602)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C][CH][C]=C[CH2](17842)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C[CH][CH][C][CH2](17843)'],
    products = ['[CH2][C]=CC=[C][CH2](16912)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['C=C1[CH][CH][C]C1(17844)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.98918e+10,'s^-1'), n=0.314867, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[C]1[CH]C=[C]CC1(17845)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.26053e+10,'s^-1'), n=0.2505, Ea=(112.349,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #4134',
    isomers = [
        '[CH2][C]=CC=[C][CH2](16912)',
    ],
    reactants = [
        ('C3H3(5450)', 'C3H3(5450)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4134',
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


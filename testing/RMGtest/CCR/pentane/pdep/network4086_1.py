species(
    label = '[CH]C(=C)C([CH2])C=[CH](16347)',
    structure = SMILES('[CH]C(=C)C([CH2])C=[CH]'),
    E0 = (838.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,408.879,408.879,408.879,408.88],'cm^-1')),
        HinderedRotor(inertia=(0.451149,'amu*angstrom^2'), symmetry=1, barrier=(53.5224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.451147,'amu*angstrom^2'), symmetry=1, barrier=(53.5224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.451143,'amu*angstrom^2'), symmetry=1, barrier=(53.5224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.451151,'amu*angstrom^2'), symmetry=1, barrier=(53.5224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676133,0.0742583,-6.93413e-05,3.84038e-08,-9.03845e-12,101026,28.3972], Tmin=(100,'K'), Tmax=(1005.16,'K')), NASAPolynomial(coeffs=[9.70687,0.0383214,-1.57135e-05,2.83603e-09,-1.92288e-13,99210.2,-15.2175], Tmin=(1005.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(838.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]C=CCC([CH])=C(16390)',
    structure = SMILES('[CH]C=CCC([CH])=C'),
    E0 = (755.038,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624932,0.064753,-3.00173e-05,7.57281e-10,2.36403e-12,90939.4,29.7277], Tmin=(100,'K'), Tmax=(1180.86,'K')), NASAPolynomial(coeffs=[11.9939,0.04009,-1.6279e-05,2.93185e-09,-1.98763e-13,87288.9,-31.099], Tmin=(1180.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)[CH]CC=[CH](18371)',
    structure = SMILES('[CH]C([CH2])=CCC=[CH]'),
    E0 = (782.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.589258,0.0654944,-3.76798e-05,4.95311e-09,2.06753e-12,94297.8,29.4251], Tmin=(100,'K'), Tmax=(1099.83,'K')), NASAPolynomial(coeffs=[13.8548,0.0322737,-1.28638e-05,2.33207e-09,-1.60122e-13,90471.1,-39.9673], Tmin=(1099.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]=C([CH2])C=C(15994)',
    structure = SMILES('[CH]C(=C)C=C'),
    E0 = (427.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.11781,'amu*angstrom^2'), symmetry=1, barrier=(48.6927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11753,'amu*angstrom^2'), symmetry=1, barrier=(48.6861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85302,0.0366828,7.83863e-06,-3.7032e-08,1.71677e-11,51444.9,17.4107], Tmin=(100,'K'), Tmax=(958.444,'K')), NASAPolynomial(coeffs=[11.8211,0.0206752,-7.16379e-06,1.26433e-09,-8.87621e-14,48358.6,-36.3891], Tmin=(958.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH2])=CC=[CH](17869)',
    structure = SMILES('[CH]C([CH2])=CC=[CH]'),
    E0 = (754.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650,196.257,196.342,196.517],'cm^-1')),
        HinderedRotor(inertia=(1.82611,'amu*angstrom^2'), symmetry=1, barrier=(50.2313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82837,'amu*angstrom^2'), symmetry=1, barrier=(50.2371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84001,'amu*angstrom^2'), symmetry=1, barrier=(50.2315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39637,0.047072,-6.79706e-06,-2.63284e-08,1.4179e-11,90869.6,22.9728], Tmin=(100,'K'), Tmax=(955.335,'K')), NASAPolynomial(coeffs=[12.8884,0.023431,-8.10907e-06,1.4059e-09,-9.68762e-14,87556.9,-37.7906], Tmin=(955.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]C(=C)C([CH])C=[CH](18372)',
    structure = SMILES('[CH]C(=C)C([CH])C=[CH]'),
    E0 = (1082.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71798,0.0742792,-7.15328e-05,3.92128e-08,-9.05109e-12,130266,27.4715], Tmin=(100,'K'), Tmax=(1024.83,'K')), NASAPolynomial(coeffs=[10.4775,0.0361871,-1.57794e-05,2.94469e-09,-2.03795e-13,128265,-19.8521], Tmin=(1024.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1082.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=CC([CH2])C([CH])=C(18373)',
    structure = SMILES('[C]=CC([CH2])C([CH])=C'),
    E0 = (1149.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05403,0.0718901,-6.42159e-05,1.56606e-08,1.19528e-11,138412,27.083], Tmin=(100,'K'), Tmax=(587.289,'K')), NASAPolynomial(coeffs=[7.93881,0.0391851,-1.69186e-05,3.10321e-09,-2.11143e-13,137358,-4.55077], Tmin=(587.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1149.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C1C=CC1(17978)',
    structure = SMILES('[CH]C(=C)C1C=CC1'),
    E0 = (547.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,350,440,435,1725,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20891,0.0485725,8.29704e-06,-4.04799e-08,1.75514e-11,65996.7,22.824], Tmin=(100,'K'), Tmax=(1010.51,'K')), NASAPolynomial(coeffs=[12.3583,0.0337604,-1.3241e-05,2.44458e-09,-1.7226e-13,62246.3,-38.4895], Tmin=(1010.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(=C)C=C(17369)',
    structure = SMILES('[CH]C(=C)C(=C)C=C'),
    E0 = (497.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,2980,3010,3040,3070,3100,1330,1380,1430,900,975,1050,1000,1025,1050,1600,1650,1700,290.345,291.77,291.885,294.194,299.961],'cm^-1')),
        HinderedRotor(inertia=(0.791377,'amu*angstrom^2'), symmetry=1, barrier=(49.8993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852343,'amu*angstrom^2'), symmetry=1, barrier=(49.9678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764697,'amu*angstrom^2'), symmetry=1, barrier=(49.9941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979797,0.0546567,-6.05821e-06,-2.81639e-08,1.39585e-11,60001.1,23.5828], Tmin=(100,'K'), Tmax=(1001.96,'K')), NASAPolynomial(coeffs=[13.0541,0.032704,-1.24918e-05,2.26433e-09,-1.57811e-13,56263.8,-41.2679], Tmin=(1001.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[C]([CH2])C=[CH](18374)',
    structure = SMILES('[CH]C=C([CH2])C([CH])[CH2]'),
    E0 = (1056.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.408045,0.0709957,-4.90398e-05,1.27408e-08,6.44577e-13,127181,31.7073], Tmin=(100,'K'), Tmax=(1013.93,'K')), NASAPolynomial(coeffs=[14.0629,0.0314066,-1.15977e-05,2.01271e-09,-1.35126e-13,123678,-37.9782], Tmin=(1013.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1056.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C([CH2])[C]=[CH](18375)',
    structure = SMILES('[CH]C([CH2])C([CH2])[C]=[CH]'),
    E0 = (1241.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.447356,0.0758514,-8.221e-05,4.9437e-08,-1.18905e-11,149483,34.4811], Tmin=(100,'K'), Tmax=(1016.23,'K')), NASAPolynomial(coeffs=[13.0276,0.0263338,-9.11921e-06,1.48754e-09,-9.45254e-14,146926,-26.4136], Tmin=(1016.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1241.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C](C)C([CH2])[C]=[CH](18376)',
    structure = SMILES('[CH][C](C)C([CH2])[C]=[CH]'),
    E0 = (1222.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,269.285,887.646,1249.94,1914.88],'cm^-1')),
        HinderedRotor(inertia=(0.095706,'amu*angstrom^2'), symmetry=1, barrier=(3.33382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095706,'amu*angstrom^2'), symmetry=1, barrier=(3.33382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095706,'amu*angstrom^2'), symmetry=1, barrier=(3.33382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095706,'amu*angstrom^2'), symmetry=1, barrier=(3.33382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095706,'amu*angstrom^2'), symmetry=1, barrier=(3.33382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640324,0.080098,-0.000111386,9.23891e-08,-3.06833e-11,147103,33.1437], Tmin=(100,'K'), Tmax=(862.982,'K')), NASAPolynomial(coeffs=[7.43373,0.0368857,-1.58973e-05,2.88007e-09,-1.92551e-13,146367,3.90001], Tmin=(862.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1222.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Cds_S) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]1CCC1C=[CH](18377)',
    structure = SMILES('[CH][C]1CCC1C=[CH]'),
    E0 = (934.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45167,0.0449569,2.01405e-06,-2.99286e-08,1.32304e-11,112440,27.5449], Tmin=(100,'K'), Tmax=(1032.86,'K')), NASAPolynomial(coeffs=[11.7248,0.0298703,-1.19451e-05,2.23333e-09,-1.58139e-13,109000,-28.7264], Tmin=(1032.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(934.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]1CC=CC1[CH2](18378)',
    structure = SMILES('[CH][C]1CC=CC1[CH2]'),
    E0 = (796.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65535,0.0406451,1.13978e-05,-3.93044e-08,1.69066e-11,95848,26.272], Tmin=(100,'K'), Tmax=(989.713,'K')), NASAPolynomial(coeffs=[10.6853,0.029799,-1.10374e-05,1.99289e-09,-1.39299e-13,92804.4,-23.5452], Tmin=(989.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(796.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Tertalkyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])CC1C=[CH](18379)',
    structure = SMILES('[CH]C1([CH2])CC1C=[CH]'),
    E0 = (945.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979162,0.0484858,1.55029e-05,-6.09436e-08,2.82461e-11,113899,27.1099], Tmin=(100,'K'), Tmax=(963.933,'K')), NASAPolynomial(coeffs=[18.75,0.018556,-6.10117e-06,1.15111e-09,-8.77488e-14,108437,-68.5294], Tmin=(963.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(945.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Neopentyl) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1([CH2])C=CC1[CH2](18380)',
    structure = SMILES('[CH]C1([CH2])C=CC1[CH2]'),
    E0 = (906.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18034,0.0444586,2.39261e-05,-6.8094e-08,3.07552e-11,109151,26.7129], Tmin=(100,'K'), Tmax=(951.377,'K')), NASAPolynomial(coeffs=[17.3983,0.0198713,-6.05024e-06,1.08179e-09,-8.0637e-14,104092,-61.0909], Tmin=(951.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(906.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Neopentyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(=C)C=[CH](18381)',
    structure = SMILES('[CH]C(=C)C(=C)C=[CH]'),
    E0 = (744.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3120,650,792.5,1650,229.22,229.221,229.244,229.308],'cm^-1')),
        HinderedRotor(inertia=(1.33851,'amu*angstrom^2'), symmetry=1, barrier=(49.9062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33807,'amu*angstrom^2'), symmetry=1, barrier=(49.9063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33828,'amu*angstrom^2'), symmetry=1, barrier=(49.9059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897875,0.0582752,-2.30465e-05,-1.03626e-08,7.86993e-12,89721.2,24.28], Tmin=(100,'K'), Tmax=(1013.69,'K')), NASAPolynomial(coeffs=[13.3429,0.0297897,-1.14112e-05,2.05464e-09,-1.42042e-13,86138.6,-41.1553], Tmin=(1013.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C([CH2])C#C(18382)',
    structure = SMILES('[CH]C(=C)C([CH2])C#C'),
    E0 = (733.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,288.211,288.214,288.214,288.215],'cm^-1')),
        HinderedRotor(inertia=(0.87112,'amu*angstrom^2'), symmetry=1, barrier=(51.3489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871113,'amu*angstrom^2'), symmetry=1, barrier=(51.3489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.871112,'amu*angstrom^2'), symmetry=1, barrier=(51.3489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87112,'amu*angstrom^2'), symmetry=1, barrier=(51.3489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507458,0.0758639,-7.59496e-05,4.30857e-08,-1.00295e-11,88339.3,25.9282], Tmin=(100,'K'), Tmax=(1034.59,'K')), NASAPolynomial(coeffs=[11.8326,0.0320782,-1.24675e-05,2.17968e-09,-1.44957e-13,85995.9,-29.0944], Tmin=(1034.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C(=C)[CH][CH2](16938)',
    structure = SMILES('[CH]C([CH2])=C[CH2]'),
    E0 = (604.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,392.587,394.013,394.853],'cm^-1')),
        HinderedRotor(inertia=(0.444068,'amu*angstrom^2'), symmetry=1, barrier=(50.6114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460465,'amu*angstrom^2'), symmetry=1, barrier=(50.6147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.463464,'amu*angstrom^2'), symmetry=1, barrier=(50.6331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81612,0.0388093,1.79373e-07,-2.59922e-08,1.22402e-11,72818.5,19.2374], Tmin=(100,'K'), Tmax=(992.532,'K')), NASAPolynomial(coeffs=[10.7608,0.0238542,-9.09623e-06,1.64932e-09,-1.15257e-13,70004,-29.0821], Tmin=(992.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(604.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C(=C)[C]([CH2])C=[CH](18383)',
    structure = SMILES('[CH]C([CH2])=C([CH2])C=[CH]'),
    E0 = (868.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,311.845,311.848,311.851],'cm^-1')),
        HinderedRotor(inertia=(0.730015,'amu*angstrom^2'), symmetry=1, barrier=(50.3788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.73001,'amu*angstrom^2'), symmetry=1, barrier=(50.3788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.730006,'amu*angstrom^2'), symmetry=1, barrier=(50.3788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.730026,'amu*angstrom^2'), symmetry=1, barrier=(50.3788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562227,0.0644614,-2.83536e-05,-1.3255e-08,1.10916e-11,104601,26.8479], Tmin=(100,'K'), Tmax=(960.401,'K')), NASAPolynomial(coeffs=[16.0075,0.0261366,-9.1099e-06,1.5791e-09,-1.08491e-13,100435,-53.2861], Tmin=(960.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(868.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C([CH2])[C]=[CH](18384)',
    structure = SMILES('[CH]C(=C)C([CH2])[C]=[CH]'),
    E0 = (1076.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,424.009,424.019,424.03,424.056],'cm^-1')),
        HinderedRotor(inertia=(0.423961,'amu*angstrom^2'), symmetry=1, barrier=(54.0979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423994,'amu*angstrom^2'), symmetry=1, barrier=(54.0975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.424089,'amu*angstrom^2'), symmetry=1, barrier=(54.0978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.424022,'amu*angstrom^2'), symmetry=1, barrier=(54.0976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.535679,0.0801753,-9.90148e-05,7.45211e-08,-2.3209e-11,129634,29.3718], Tmin=(100,'K'), Tmax=(840.707,'K')), NASAPolynomial(coeffs=[8.62107,0.0376463,-1.58909e-05,2.86156e-09,-1.91682e-13,128418,-7.37939], Tmin=(840.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1076.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C([CH2])C=[CH](18385)',
    structure = SMILES('[CH]C(=[CH])C([CH2])C=[CH]'),
    E0 = (1086.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3115,3125,620,680,785,800,1600,1700,367.689,367.689,367.689],'cm^-1')),
        HinderedRotor(inertia=(0.556103,'amu*angstrom^2'), symmetry=1, barrier=(53.3511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.556103,'amu*angstrom^2'), symmetry=1, barrier=(53.3511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.556103,'amu*angstrom^2'), symmetry=1, barrier=(53.3511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.556104,'amu*angstrom^2'), symmetry=1, barrier=(53.3511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550299,0.0783807,-8.80417e-05,5.8397e-08,-1.60701e-11,130748,28.5597], Tmin=(100,'K'), Tmax=(877.152,'K')), NASAPolynomial(coeffs=[9.95646,0.0354847,-1.46827e-05,2.63903e-09,-1.77631e-13,129098,-15.5862], Tmin=(877.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1086.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)[C](C)C=[CH](18386)',
    structure = SMILES('[CH]C([CH2])=C(C)C=[CH]'),
    E0 = (717.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518899,0.0674717,-3.87246e-05,2.48734e-09,4.05323e-12,86378.8,26.7228], Tmin=(100,'K'), Tmax=(1009.93,'K')), NASAPolynomial(coeffs=[14.029,0.0316606,-1.18224e-05,2.08072e-09,-1.41406e-13,82747.4,-43.0578], Tmin=(1009.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C([CH2])[C]=C(17371)',
    structure = SMILES('[CH]C(=C)C([CH2])[C]=C'),
    E0 = (829.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716517,0.0752951,-7.71454e-05,4.96723e-08,-1.37709e-11,99909.4,28.3259], Tmin=(100,'K'), Tmax=(858.037,'K')), NASAPolynomial(coeffs=[8.06558,0.0410356,-1.72546e-05,3.13972e-09,-2.13231e-13,98648.2,-6.00392], Tmin=(858.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(829.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(C)[C]=[CH](18387)',
    structure = SMILES('[CH]C(=C)C(C)[C]=[CH]'),
    E0 = (871.747,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,319.789,319.791,319.791,319.792],'cm^-1')),
        HinderedRotor(inertia=(0.722681,'amu*angstrom^2'), symmetry=1, barrier=(52.4441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722669,'amu*angstrom^2'), symmetry=1, barrier=(52.4441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722672,'amu*angstrom^2'), symmetry=1, barrier=(52.4441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722665,'amu*angstrom^2'), symmetry=1, barrier=(52.4441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759694,0.0754667,-7.51522e-05,4.61801e-08,-1.23524e-11,104960,26.798], Tmin=(100,'K'), Tmax=(880.507,'K')), NASAPolynomial(coeffs=[8.03729,0.0424061,-1.8832e-05,3.53836e-09,-2.45417e-13,103678,-7.3862], Tmin=(880.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(871.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)[C]([CH2])C=C(17370)',
    structure = SMILES('[CH]C([CH2])=C([CH2])C=C'),
    E0 = (621.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,358.386,358.386,358.387,358.387],'cm^-1')),
        HinderedRotor(inertia=(0.553787,'amu*angstrom^2'), symmetry=1, barrier=(50.4747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553788,'amu*angstrom^2'), symmetry=1, barrier=(50.4747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553788,'amu*angstrom^2'), symmetry=1, barrier=(50.4747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.553788,'amu*angstrom^2'), symmetry=1, barrier=(50.4747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634741,0.0609528,-1.17418e-05,-3.05859e-08,1.69897e-11,74880.8,26.1845], Tmin=(100,'K'), Tmax=(961.787,'K')), NASAPolynomial(coeffs=[15.7722,0.0289614,-1.01394e-05,1.77681e-09,-1.23273e-13,70536.9,-53.7003], Tmin=(961.787,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C(C)C=[CH](18388)',
    structure = SMILES('[CH]C(=[CH])C(C)C=[CH]'),
    E0 = (881.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,236.305,236.305,236.305],'cm^-1')),
        HinderedRotor(inertia=(1.29863,'amu*angstrom^2'), symmetry=1, barrier=(51.4581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29862,'amu*angstrom^2'), symmetry=1, barrier=(51.4581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29862,'amu*angstrom^2'), symmetry=1, barrier=(51.4581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29862,'amu*angstrom^2'), symmetry=1, barrier=(51.4581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.720979,0.0744594,-6.76648e-05,3.55856e-08,-8.00546e-12,106076,26.1666], Tmin=(100,'K'), Tmax=(1040.34,'K')), NASAPolynomial(coeffs=[9.83881,0.0394013,-1.71154e-05,3.1918e-09,-2.20825e-13,104179,-18.1819], Tmin=(1040.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(881.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C([CH2])C=C(17372)',
    structure = SMILES('[CH]C(=[CH])C([CH2])C=C'),
    E0 = (838.988,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,408.878,408.879,408.88,408.881],'cm^-1')),
        HinderedRotor(inertia=(0.451144,'amu*angstrom^2'), symmetry=1, barrier=(53.5224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45115,'amu*angstrom^2'), symmetry=1, barrier=(53.5223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.451149,'amu*angstrom^2'), symmetry=1, barrier=(53.5225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45115,'amu*angstrom^2'), symmetry=1, barrier=(53.5224,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676133,0.0742583,-6.93413e-05,3.84038e-08,-9.03845e-12,101026,27.704], Tmin=(100,'K'), Tmax=(1005.16,'K')), NASAPolynomial(coeffs=[9.70687,0.0383214,-1.57135e-05,2.83603e-09,-1.92288e-13,99210.2,-15.9107], Tmin=(1005.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(838.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3375.81,'J/mol'), sigma=(5.99702,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.29 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778217,0.0684668,-6.41125e-05,3.31395e-08,-6.99132e-12,116670,31.6665], Tmin=(100,'K'), Tmax=(1137.47,'K')), NASAPolynomial(coeffs=[12.2659,0.0280694,-1.08396e-05,1.91637e-09,-1.28891e-13,114056,-25.2345], Tmin=(1137.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(969.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH2])C=[CH](14965)',
    structure = SMILES('[CH]=[C]C([CH2])C=[CH]'),
    E0 = (992.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(0.346294,'amu*angstrom^2'), symmetry=1, barrier=(7.96199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346556,'amu*angstrom^2'), symmetry=1, barrier=(7.96801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.877943,'amu*angstrom^2'), symmetry=1, barrier=(20.1856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22561,0.0647448,-8.904e-05,7.01424e-08,-2.20432e-11,119442,25.4311], Tmin=(100,'K'), Tmax=(875.088,'K')), NASAPolynomial(coeffs=[8.41707,0.0254017,-1.05092e-05,1.8651e-09,-1.23173e-13,118431,-6.88833], Tmin=(875.088,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(992.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[C]=C([CH2])C([CH2])C=[CH](18389)',
    structure = SMILES('[C]C(=C)C([CH2])C=[CH]'),
    E0 = (1137.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,238.981,239.873],'cm^-1')),
        HinderedRotor(inertia=(0.246827,'amu*angstrom^2'), symmetry=1, barrier=(10.0055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693604,'amu*angstrom^2'), symmetry=1, barrier=(28.2422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245493,'amu*angstrom^2'), symmetry=1, barrier=(10.005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562595,0.0771706,-9.23864e-05,6.06267e-08,-1.59714e-11,136966,26.1487], Tmin=(100,'K'), Tmax=(926.068,'K')), NASAPolynomial(coeffs=[12.2136,0.0268472,-1.08768e-05,1.94998e-09,-1.31419e-13,134808,-29.1659], Tmin=(926.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1137.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P) + radical(CJ3)"""),
)

species(
    label = '[CH]=CC([CH2])C1=CC1(18390)',
    structure = SMILES('[CH]=CC([CH2])C1=CC1'),
    E0 = (736.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.608342,0.0675963,-5.90598e-05,2.77064e-08,-5.23099e-12,88748.3,27.1159], Tmin=(100,'K'), Tmax=(1276.09,'K')), NASAPolynomial(coeffs=[14.2174,0.0249383,-8.91752e-06,1.51092e-09,-9.90849e-14,85275,-41.8579], Tmin=(1276.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(736.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1CCC1=[CH](18391)',
    structure = SMILES('[CH]=CC1CCC1=[CH]'),
    E0 = (681.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11281,0.0550763,-2.9018e-05,3.4616e-09,1.23907e-12,82125.5,23.7781], Tmin=(100,'K'), Tmax=(1216.11,'K')), NASAPolynomial(coeffs=[12.4304,0.0296652,-1.22471e-05,2.2562e-09,-1.55314e-13,78499.2,-36.6288], Tmin=(1216.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1CC=C1[CH2](18349)',
    structure = SMILES('[CH]=CC1C[CH]C1=C'),
    E0 = (575.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45158,0.0456444,-1.51811e-06,-2.41025e-08,1.05173e-11,69368.3,21.0458], Tmin=(100,'K'), Tmax=(1074.4,'K')), NASAPolynomial(coeffs=[11.2271,0.0316098,-1.31409e-05,2.47942e-09,-1.75187e-13,65977.2,-32.8228], Tmin=(1074.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC=CC1[CH2](18392)',
    structure = SMILES('[CH]=C1CC=CC1[CH2]'),
    E0 = (533.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60763,0.0421388,6.81368e-06,-3.39405e-08,1.4776e-11,64318.2,23.6684], Tmin=(100,'K'), Tmax=(1003.8,'K')), NASAPolynomial(coeffs=[10.6215,0.0302919,-1.14546e-05,2.08239e-09,-1.45569e-13,61295.8,-25.8933], Tmin=(1003.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=CC=CC1[CH2](18393)',
    structure = SMILES('[CH2]C1=CC=CC1[CH2]'),
    E0 = (372.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54132,0.0390304,2.98789e-05,-6.75948e-08,2.93662e-11,44932.9,20.9616], Tmin=(100,'K'), Tmax=(943.862,'K')), NASAPolynomial(coeffs=[13.7154,0.0251091,-7.86484e-06,1.34992e-09,-9.53415e-14,40956.8,-45.9569], Tmin=(943.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclopentadiene) + radical(Isobutyl) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC(=C)C(=[CH])C(18394)',
    structure = SMILES('[CH]=CC(=C)C(=[CH])C'),
    E0 = (621.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.783384,0.0623427,-4.16789e-05,9.58346e-09,8.79666e-13,74859.6,23.9862], Tmin=(100,'K'), Tmax=(1063.94,'K')), NASAPolynomial(coeffs=[14.0143,0.0262438,-1.00208e-05,1.79996e-09,-1.23699e-13,71272,-44.2943], Tmin=(1063.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(=C)C([CH2])=C(16358)',
    structure = SMILES('[CH]=CC(=C)C([CH2])=C'),
    E0 = (525.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3120,650,792.5,1650,312.221],'cm^-1')),
        HinderedRotor(inertia=(0.419541,'amu*angstrom^2'), symmetry=1, barrier=(29.022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419531,'amu*angstrom^2'), symmetry=1, barrier=(29.0215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419599,'amu*angstrom^2'), symmetry=1, barrier=(29.0222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916807,0.0556424,-1.41974e-05,-2.39115e-08,1.38858e-11,63360.8,22.6901], Tmin=(100,'K'), Tmax=(978.291,'K')), NASAPolynomial(coeffs=[15.5423,0.0238956,-8.53424e-06,1.54144e-09,-1.09266e-13,59156.8,-54.4094], Tmin=(978.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C#C(14975)',
    structure = SMILES('[CH]=CC([CH2])C#C'),
    E0 = (648.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2175,525,3120,650,792.5,1650,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(2.46263,'amu*angstrom^2'), symmetry=1, barrier=(56.6206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843691,'amu*angstrom^2'), symmetry=1, barrier=(19.3981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843625,'amu*angstrom^2'), symmetry=1, barrier=(19.3966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16694,0.060797,-6.72891e-05,4.05297e-08,-9.72226e-12,78149.1,22.0966], Tmin=(100,'K'), Tmax=(1019.62,'K')), NASAPolynomial(coeffs=[11.5905,0.0199042,-7.12909e-06,1.19394e-09,-7.73697e-14,76023.5,-28.3934], Tmin=(1019.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]([CH2])C(=[CH])C(18395)',
    structure = SMILES('[CH]C(C)=C([CH2])C=[CH]'),
    E0 = (750.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.212858,0.0741669,-5.43381e-05,1.61029e-08,-3.13713e-13,90412.1,26.5447], Tmin=(100,'K'), Tmax=(1046.86,'K')), NASAPolynomial(coeffs=[15.7744,0.0293665,-1.11507e-05,1.97647e-09,-1.34563e-13,86350.6,-53.0802], Tmin=(1046.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(750.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]([CH2])C([CH2])=C(16360)',
    structure = SMILES('[CH]=CC([CH2])=C([CH2])[CH2]'),
    E0 = (615.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.235738,'amu*angstrom^2'), symmetry=1, barrier=(53.5704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235795,'amu*angstrom^2'), symmetry=1, barrier=(53.5691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235637,'amu*angstrom^2'), symmetry=1, barrier=(53.5716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235586,'amu*angstrom^2'), symmetry=1, barrier=(53.5698,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847894,0.0555641,-5.21339e-06,-3.90291e-08,2.10496e-11,74208.6,26.2721], Tmin=(100,'K'), Tmax=(938.152,'K')), NASAPolynomial(coeffs=[16.8129,0.0219552,-6.57589e-06,1.09362e-09,-7.62809e-14,69696.5,-57.8131], Tmin=(938.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH2])C(=[CH])C(18396)',
    structure = SMILES('[CH]=[C]C([CH2])C(=[CH])C'),
    E0 = (953.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(0.33547,'amu*angstrom^2'), symmetry=1, barrier=(7.71312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335513,'amu*angstrom^2'), symmetry=1, barrier=(7.7141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.335411,'amu*angstrom^2'), symmetry=1, barrier=(7.71176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757396,'amu*angstrom^2'), symmetry=1, barrier=(17.414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.422753,0.0842677,-0.000117955,9.52132e-08,-3.06825e-11,114772,29.0697], Tmin=(100,'K'), Tmax=(864.019,'K')), NASAPolynomial(coeffs=[9.13031,0.0343566,-1.46403e-05,2.63852e-09,-1.75877e-13,113626,-9.59309], Tmin=(864.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(953.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH2])C([CH2])=C(16361)',
    structure = SMILES('[CH]=[C]C([CH2])C([CH2])=C'),
    E0 = (857.644,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,223.023],'cm^-1')),
        HinderedRotor(inertia=(2.05651,'amu*angstrom^2'), symmetry=1, barrier=(74.5446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00236772,'amu*angstrom^2'), symmetry=1, barrier=(13.7264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387071,'amu*angstrom^2'), symmetry=1, barrier=(13.7126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08996,'amu*angstrom^2'), symmetry=1, barrier=(74.5414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.597208,0.0769675,-8.77711e-05,5.7218e-08,-1.52553e-11,103271,27.6335], Tmin=(100,'K'), Tmax=(908.803,'K')), NASAPolynomial(coeffs=[10.8963,0.0316375,-1.29537e-05,2.33516e-09,-1.57874e-13,101399,-21.069], Tmin=(908.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(857.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (838.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (998.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (998.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (838.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1192.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1293.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1361.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (847.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (902.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1079.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1305.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1247.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (964.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (877.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (945.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (906.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (974.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (965.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (974.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1039.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (850.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1108.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1195.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1080.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1288.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1297.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (937.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (944.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1023.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1031.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (925.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (892.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (1214.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1430.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1349.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (847.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (847.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (847.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (846.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (846.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (902.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (902.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (1047.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (862.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (1041.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (1031.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (1055.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (1251.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['C3H3(5450)', 'CH2CHCHCH(4849)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C=CCC([CH])=C(16390)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C(=C)[CH]CC=[CH](18371)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['C2H2(1342)', '[CH]=C([CH2])C=C(15994)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH]C([CH2])=CC=[CH](17869)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C(=C)C([CH])C=[CH](18372)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[C]=CC([CH2])C([CH])=C(18373)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C(=C)C1C=CC1(17978)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C(=C)C(=C)C=C(17369)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]C([CH2])[C]([CH2])C=[CH](18374)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]C([CH2])C([CH2])[C]=[CH](18375)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH][C](C)C([CH2])[C]=[CH](18376)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH][C]1CCC1C=[CH](18377)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH][C]1CC=CC1[CH2](18378)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.29e+09,'s^-1'), n=0.62, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C1([CH2])CC1C=[CH](18379)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(106.982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 103.0 to 107.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C1([CH2])C=CC1[CH2](18380)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(67.565,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 63.1 to 67.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH]C(=C)C(=C)C=[CH](18381)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-TwoDe_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH]C(=C)C([CH2])C#C(18382)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.276e+10,'cm^3/(mol*s)'), n=1.103, Ea=(20.5476,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 139 used for Ct-Cs_Ct-H;HJ
Exact match found for rate rule [Ct-Cs_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][CH2](16918)', 'CH2CHCHCH(4849)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00509465,'m^3/(mol*s)'), n=2.41843, Ea=(13.2728,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H2(T)(1343)', '[CH]=C([CH2])C=C(15994)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C2H2(1342)', '[CH]C(=C)[CH][CH2](16938)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C][CH2](16918)', '[CH]C=C[CH2](15066)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C2H2(T)(1343)', '[CH]C(=C)[CH][CH2](16938)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH]C(=C)[C]([CH2])C=[CH](18383)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]C(=C)C([CH2])[C]=[CH](18384)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.17499e+07,'m^3/(mol*s)'), n=0.115342, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH]C(=[CH])C([CH2])C=[CH](18385)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.63133e+13,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C(=C)[C](C)C=[CH](18386)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C(=C)C([CH2])[C]=C(17371)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(=C)C(C)[C]=[CH](18387)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]C(=C)[C]([CH2])C=C(17370)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C(=[CH])C(C)C=[CH](18388)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(222600,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]C(=[CH])C([CH2])C=C(17372)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(769414,'s^-1'), n=1.8337, Ea=(53.4313,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;XH_out] + [R5H_RSSR;Y_rad_out;Cd_H_out_singleH] for rate rule [R5H_DSSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 3.605551275463989
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C]CC([CH2])C=[CH](16339)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH2(T)(20)', '[CH]=[C]C([CH2])C=[CH](14965)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[C]=C([CH2])C([CH2])C=[CH](18389)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=CC([CH2])C1=CC1(18390)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=CC1CCC1=[CH](18391)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=CC1CC=C1[CH2](18349)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=C1CC=CC1[CH2](18392)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH2]C1=CC=CC1[CH2](18393)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_DSSD;CdsingleH_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=CC(=C)C(=[CH])C(18394)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=CC(=C)C([CH2])=C(16358)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CH2(T)(20)', '[CH]=CC([CH2])C#C(14975)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C3H3(5450)', '[CH]C=C[CH2](15066)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(0.00505236,'m^3/(mol*s)'), n=2.48091, Ea=(29.9954,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cdd_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=C[C]([CH2])C(=[CH])C(18395)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    products = ['[CH]=C[C]([CH2])C([CH2])=C(16360)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C]C([CH2])C(=[CH])C(18396)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.35985e+07,'s^-1'), n=1.33667, Ea=(102.123,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=[C]C([CH2])C([CH2])=C(16361)'],
    products = ['[CH]C(=C)C([CH2])C=[CH](16347)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(3.60477e+12,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4086',
    isomers = [
        '[CH]C(=C)C([CH2])C=[CH](16347)',
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
    label = 'PDepNetwork #4086',
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


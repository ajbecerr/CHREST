species(
    label = '[CH]=[C]CC(=[CH])[CH2](16901)',
    structure = SMILES('[CH]C(=C)C[C]=[CH]'),
    E0 = (905.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,351.264,352.028,352.055,352.809],'cm^-1')),
        HinderedRotor(inertia=(0.584914,'amu*angstrom^2'), symmetry=1, barrier=(51.0487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585247,'amu*angstrom^2'), symmetry=1, barrier=(51.0279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58312,'amu*angstrom^2'), symmetry=1, barrier=(51.0421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46847,0.0539598,-4.00949e-05,1.60924e-08,-2.72863e-12,108977,24.6512], Tmin=(100,'K'), Tmax=(1348.6,'K')), NASAPolynomial(coeffs=[10.0074,0.028633,-1.19247e-05,2.16667e-09,-1.47114e-13,106674,-19.0978], Tmin=(1348.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]=C([CH2])C(=[CH])[CH2](16900)',
    structure = SMILES('[CH]C(=C)C([CH])=C'),
    E0 = (760.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00411,0.0566641,-2.24676e-05,-1.04642e-08,8.05012e-12,91536,21.474], Tmin=(100,'K'), Tmax=(992.534,'K')), NASAPolynomial(coeffs=[12.6805,0.0293738,-1.10972e-05,1.96361e-09,-1.34205e-13,88244.5,-39.6749], Tmin=(992.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=[CH](9646)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([554.803,1738.79,3454.47],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97621,0.00212915,-8.08978e-08,-3.83305e-10,9.76908e-14,101881,6.00119], Tmin=(100,'K'), Tmax=(1982.31,'K')), NASAPolynomial(coeffs=[5.05695,0.00131032,-4.91873e-07,1.01502e-10,-7.16167e-15,101185,-0.627232], Tmin=(1982.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])[CH2](16419)',
    structure = SMILES('[CH]C([CH2])=C'),
    E0 = (488.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,412.019,412.344,412.414],'cm^-1')),
        HinderedRotor(inertia=(0.422651,'amu*angstrom^2'), symmetry=1, barrier=(50.9185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421927,'amu*angstrom^2'), symmetry=1, barrier=(50.9131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54708,0.0237852,1.61126e-05,-3.4475e-08,1.40391e-11,58778.2,14.8721], Tmin=(100,'K'), Tmax=(988.925,'K')), NASAPolynomial(coeffs=[8.49332,0.0194138,-7.10753e-06,1.30164e-09,-9.24087e-14,56639.8,-18.6144], Tmin=(988.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    label = '[C]=[C]CC([CH])=C(17889)',
    structure = SMILES('[C]=[C]CC([CH])=C'),
    E0 = (1216.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,406.475,407.457,408.747,409.513,411.574],'cm^-1')),
        HinderedRotor(inertia=(0.431261,'amu*angstrom^2'), symmetry=1, barrier=(51.8595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447797,'amu*angstrom^2'), symmetry=1, barrier=(51.9518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.442958,'amu*angstrom^2'), symmetry=1, barrier=(51.9761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66407,0.054766,-5.27754e-05,3.1249e-08,-8.14771e-12,146371,23.938], Tmin=(100,'K'), Tmax=(895.476,'K')), NASAPolynomial(coeffs=[6.77979,0.0319148,-1.44981e-05,2.75248e-09,-1.92113e-13,145455,-0.177699], Tmin=(895.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1216.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[CH][C]=[CH](17890)',
    structure = SMILES('[CH][C]=CC([CH])[CH2]'),
    E0 = (1177.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32031,0.0568939,-4.76659e-05,2.29218e-08,-4.64137e-12,141722,27.8696], Tmin=(100,'K'), Tmax=(1162.84,'K')), NASAPolynomial(coeffs=[9.57162,0.0285101,-1.10517e-05,1.93013e-09,-1.28274e-13,139803,-13.1828], Tmin=(1162.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1177.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1CC(=[CH])C1(17891)',
    structure = SMILES('[CH][C]1CC(=[CH])C1'),
    E0 = (963.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12622,0.0397992,-2.11935e-05,5.13204e-09,-4.84046e-13,116002,21.5153], Tmin=(100,'K'), Tmax=(2415.38,'K')), NASAPolynomial(coeffs=[16.1781,0.016529,-6.74251e-06,1.14352e-09,-7.12292e-14,109214,-58.6689], Tmin=(2415.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(963.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Tertalkyl) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]1C[C]=CC1(17892)',
    structure = SMILES('[CH][C]1C[C]=CC1'),
    E0 = (880.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26552,0.0311673,3.26431e-06,-1.91547e-08,7.51345e-12,105946,21.4729], Tmin=(100,'K'), Tmax=(1104.61,'K')), NASAPolynomial(coeffs=[7.88246,0.0262325,-1.09538e-05,2.052e-09,-1.4362e-13,103766,-10.4386], Tmin=(1104.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(880.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Tertalkyl) + radical(cyclopentene-vinyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])CC1=[CH](17893)',
    structure = SMILES('[CH]C1([CH2])CC1=[CH]'),
    E0 = (1022.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40888,0.0447839,-7.06417e-06,-2.78811e-08,1.50218e-11,123118,22.0841], Tmin=(100,'K'), Tmax=(976.618,'K')), NASAPolynomial(coeffs=[15.7038,0.0149349,-5.2985e-06,1.00371e-09,-7.49006e-14,118957,-53.5493], Tmin=(976.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1022.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Neopentyl) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1([CH2])C=[C]C1(17875)',
    structure = SMILES('[CH]C1([CH2])C=[C]C1'),
    E0 = (985.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80179,0.0349229,1.56066e-05,-4.7075e-08,2.06882e-11,118644,21.8688], Tmin=(100,'K'), Tmax=(982.923,'K')), NASAPolynomial(coeffs=[14.1129,0.0170824,-6.39784e-06,1.23981e-09,-9.29801e-14,114666,-45.2405], Tmin=(982.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(985.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Neopentyl) + radical(cyclobutene-vinyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])=CC#C(17865)',
    structure = SMILES('[CH]C([CH2])=CC#C'),
    E0 = (687.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,184.031,184.031,184.031],'cm^-1')),
        HinderedRotor(inertia=(2.08763,'amu*angstrom^2'), symmetry=1, barrier=(50.1717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08762,'amu*angstrom^2'), symmetry=1, barrier=(50.1717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08762,'amu*angstrom^2'), symmetry=1, barrier=(50.1717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40126,0.0510272,-3.38163e-05,8.73285e-09,1.75219e-13,82813.3,19.5391], Tmin=(100,'K'), Tmax=(1061.46,'K')), NASAPolynomial(coeffs=[10.9855,0.0242798,-9.25892e-06,1.62518e-09,-1.09359e-13,80250.8,-29.7573], Tmin=(1061.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH2])=C[C]=[CH](17867)',
    structure = SMILES('[CH]C([CH2])=C[C]=[CH]'),
    E0 = (953.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.20061,'amu*angstrom^2'), symmetry=1, barrier=(50.5962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20252,'amu*angstrom^2'), symmetry=1, barrier=(50.6403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20003,'amu*angstrom^2'), symmetry=1, barrier=(50.5831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33019,0.0534685,-3.53423e-05,5.78722e-09,2.75888e-12,114801,22.8075], Tmin=(100,'K'), Tmax=(930.983,'K')), NASAPolynomial(coeffs=[11.325,0.0235591,-8.15178e-06,1.354e-09,-8.87473e-14,112375,-27.7303], Tmin=(930.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(953.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C[C]=[CH](17894)',
    structure = SMILES('[CH]C(=[CH])C[C]=[CH]'),
    E0 = (1152.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,3115,3125,620,680,785,800,1600,1700,290.985,290.986,290.986],'cm^-1')),
        HinderedRotor(inertia=(0.849229,'amu*angstrom^2'), symmetry=1, barrier=(51.0264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.849227,'amu*angstrom^2'), symmetry=1, barrier=(51.0264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84923,'amu*angstrom^2'), symmetry=1, barrier=(51.0264,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49473,0.0562333,-5.21705e-05,2.74231e-08,-6.08159e-12,138693,24.2733], Tmin=(100,'K'), Tmax=(1062.37,'K')), NASAPolynomial(coeffs=[9.045,0.0278053,-1.20318e-05,2.23496e-09,-1.54248e-13,137088,-12.6092], Tmin=(1062.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1152.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[CH]=C([CH2])C=[C][CH2](16906)',
    structure = SMILES('[CH]C([CH2])=C[C]=C'),
    E0 = (706.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,218.576,218.577,218.577,218.577],'cm^-1')),
        HinderedRotor(inertia=(1.49931,'amu*angstrom^2'), symmetry=1, barrier=(50.831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49932,'amu*angstrom^2'), symmetry=1, barrier=(50.831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49932,'amu*angstrom^2'), symmetry=1, barrier=(50.831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.87,'J/mol'), sigma=(5.74688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.46 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39391,0.0500702,-1.91473e-05,-1.0955e-08,8.38317e-12,85081.9,22.1753], Tmin=(100,'K'), Tmax=(944.297,'K')), NASAPolynomial(coeffs=[11.1102,0.0263477,-9.15985e-06,1.54656e-09,-1.03096e-13,82469.5,-28.2591], Tmin=(944.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])CC=[CH](17895)',
    structure = SMILES('[CH]C(=[CH])CC=[CH]'),
    E0 = (914.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,269.115,269.132,269.21],'cm^-1')),
        HinderedRotor(inertia=(0.98311,'amu*angstrom^2'), symmetry=1, barrier=(50.4366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983977,'amu*angstrom^2'), symmetry=1, barrier=(50.439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985795,'amu*angstrom^2'), symmetry=1, barrier=(50.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16172,0.0558966,-4.19666e-05,1.64803e-08,-2.64571e-12,110105,24.9967], Tmin=(100,'K'), Tmax=(1460.72,'K')), NASAPolynomial(coeffs=[12.7024,0.0242938,-9.51404e-06,1.66906e-09,-1.10795e-13,106734,-35.0538], Tmin=(1460.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(914.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C[C]=C(17482)',
    structure = SMILES('[CH]C(=[CH])C[C]=C'),
    E0 = (905.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,351.264,352.028,352.055,352.809],'cm^-1')),
        HinderedRotor(inertia=(0.584914,'amu*angstrom^2'), symmetry=1, barrier=(51.0487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585247,'amu*angstrom^2'), symmetry=1, barrier=(51.0279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58312,'amu*angstrom^2'), symmetry=1, barrier=(51.0421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46847,0.0539598,-4.00949e-05,1.60924e-08,-2.72863e-12,108977,23.9581], Tmin=(100,'K'), Tmax=(1348.6,'K')), NASAPolynomial(coeffs=[10.0074,0.028633,-1.19247e-05,2.16667e-09,-1.47114e-13,106674,-19.7909], Tmin=(1348.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC[C]=[CH](16902)',
    structure = SMILES('[CH]=[C]CC[C]=[CH]'),
    E0 = (1033.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1700,300,440,3115,3125,620,680,785,800,1600,1700,289.879],'cm^-1')),
        HinderedRotor(inertia=(0.172772,'amu*angstrom^2'), symmetry=1, barrier=(10.4197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174352,'amu*angstrom^2'), symmetry=1, barrier=(10.4001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173831,'amu*angstrom^2'), symmetry=1, barrier=(10.3977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3204.41,'J/mol'), sigma=(5.64764,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.52 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57753,0.0567016,-6.40839e-05,4.27021e-08,-1.19594e-11,124395,24.8076], Tmin=(100,'K'), Tmax=(855.473,'K')), NASAPolynomial(coeffs=[7.8951,0.0271618,-1.22877e-05,2.33716e-09,-1.6321e-13,123314,-4.68485], Tmin=(855.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1033.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[CH]=[C]C[C]=[CH](17536)',
    structure = SMILES('[CH]=[C]C[C]=[CH]'),
    E0 = (1058.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1670,1700,300,440,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(0.550742,'amu*angstrom^2'), symmetry=1, barrier=(12.6626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551168,'amu*angstrom^2'), symmetry=1, barrier=(12.6724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19038,0.04225,-5.14744e-05,3.63597e-08,-1.06192e-11,127386,20.3858], Tmin=(100,'K'), Tmax=(827.624,'K')), NASAPolynomial(coeffs=[7.12238,0.0184138,-8.2746e-06,1.56251e-09,-1.08408e-13,126570,-2.47527], Tmin=(827.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1058.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[C]=C([CH2])C[C]=[CH](17896)',
    structure = SMILES('[C]C(=C)C[C]=[CH]'),
    E0 = (1204.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.695806,'amu*angstrom^2'), symmetry=1, barrier=(15.9979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.693712,'amu*angstrom^2'), symmetry=1, barrier=(15.9498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45544,0.0556379,-5.86935e-05,3.25234e-08,-7.23187e-12,144913,22.0468], Tmin=(100,'K'), Tmax=(1087.75,'K')), NASAPolynomial(coeffs=[11.4319,0.0189513,-8.103e-06,1.51724e-09,-1.05682e-13,142743,-26.9232], Tmin=(1087.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1204.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CJ3)"""),
)

species(
    label = '[CH]=[C]CC1=CC1(17897)',
    structure = SMILES('[CH]=[C]CC1=CC1'),
    E0 = (803.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59412,0.0449577,-2.15115e-05,-5.25818e-09,5.47779e-12,96691.4,22.6819], Tmin=(100,'K'), Tmax=(1014.19,'K')), NASAPolynomial(coeffs=[12.4287,0.0187183,-7.09631e-06,1.30094e-09,-9.16888e-14,93645.6,-33.9232], Tmin=(1014.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC(=[CH])C1(17898)',
    structure = SMILES('[CH]=C1CC(=[CH])C1'),
    E0 = (690.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19548,0.0213273,5.7517e-05,-9.36081e-08,3.82537e-11,83109,19.6133], Tmin=(100,'K'), Tmax=(950.253,'K')), NASAPolynomial(coeffs=[14.9725,0.0138637,-3.81886e-06,7.19829e-10,-5.84188e-14,78589.4,-52.3803], Tmin=(950.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1=CC(=C)C1(17860)',
    structure = SMILES('[CH]=C1[CH]C(=C)C1'),
    E0 = (544.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58998,0.00539758,0.000111832,-1.55683e-07,6.19517e-11,65617.6,18.8508], Tmin=(100,'K'), Tmax=(932.564,'K')), NASAPolynomial(coeffs=[16.8768,0.00969453,-5.57311e-07,6.51522e-11,-1.51414e-14,60101.4,-64.366], Tmin=(932.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC(=[CH])C(17899)',
    structure = SMILES('[CH]C(C)=CC#C'),
    E0 = (566.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24424,0.0530585,-3.20005e-05,6.24464e-09,7.41692e-13,68217,20.4564], Tmin=(100,'K'), Tmax=(1133.43,'K')), NASAPolynomial(coeffs=[11.8827,0.0261697,-1.05169e-05,1.90259e-09,-1.30039e-13,65121,-35.2198], Tmin=(1133.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C=CC([CH2])=C(17497)',
    structure = SMILES('C#CC=C([CH2])[CH2]'),
    E0 = (438.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(3.10171,'amu*angstrom^2'), symmetry=1, barrier=(71.3145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.10166,'amu*angstrom^2'), symmetry=1, barrier=(71.3133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142724,'amu*angstrom^2'), symmetry=1, barrier=(71.3188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39256,0.0509418,-4.21539e-05,1.91619e-08,-3.53283e-12,52833.8,18.1133], Tmin=(100,'K'), Tmax=(1303.15,'K')), NASAPolynomial(coeffs=[11.377,0.0202952,-6.87877e-06,1.11622e-09,-7.09586e-14,50231.5,-32.6999], Tmin=(1303.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(CTCC=CCJ)"""),
)

species(
    label = '[CH]=[C]CC#C(17535)',
    structure = SMILES('[CH]=[C]CC#C'),
    E0 = (743.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2175,525,3120,650,792.5,1650,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.931089,'amu*angstrom^2'), symmetry=1, barrier=(21.4076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.933124,'amu*angstrom^2'), symmetry=1, barrier=(21.4544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29144,0.037637,-3.62829e-05,1.95364e-08,-4.34211e-12,89431.9,17.0016], Tmin=(100,'K'), Tmax=(1074.29,'K')), NASAPolynomial(coeffs=[7.96109,0.0165265,-6.80666e-06,1.24437e-09,-8.52955e-14,88213.7,-10.7574], Tmin=(1074.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(C)=C[C]=[CH](17871)',
    structure = SMILES('[CH][C](C)C=C=[CH]'),
    E0 = (826.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,409.895,410.111,1935.54],'cm^-1')),
        HinderedRotor(inertia=(0.110726,'amu*angstrom^2'), symmetry=1, barrier=(13.2095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273426,'amu*angstrom^2'), symmetry=1, barrier=(72.7076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610449,'amu*angstrom^2'), symmetry=1, barrier=(72.7062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59499,0.0467614,-2.52236e-05,-2.84062e-09,5.55694e-12,99554.6,21.453], Tmin=(100,'K'), Tmax=(937.483,'K')), NASAPolynomial(coeffs=[11.3963,0.0200354,-6.61128e-06,1.09738e-09,-7.29439e-14,97053.6,-28.7372], Tmin=(937.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(826.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_T) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C][CH]C([CH2])=C(17498)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH2]'),
    E0 = (701.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.626234,'amu*angstrom^2'), symmetry=1, barrier=(53.6781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.620032,'amu*angstrom^2'), symmetry=1, barrier=(53.6613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62598,'amu*angstrom^2'), symmetry=1, barrier=(53.6644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61259,0.0445986,-1.22336e-05,-2.00637e-08,1.28141e-11,84409.4,22.2441], Tmin=(100,'K'), Tmax=(907.915,'K')), NASAPolynomial(coeffs=[12.1908,0.0192769,-5.56052e-06,8.55145e-10,-5.54365e-14,81611.4,-32.5984], Tmin=(907.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
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
    E0 = (905.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1150.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1390.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1428.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1200.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1031.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (937.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1022.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (985.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (911.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (980.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1228.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1165.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1364.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1022.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (949.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (958.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1096.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1279.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1496.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1415.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (913.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (913.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (913.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (968.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (968.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1141.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1107.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1097.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['C3H3(5450)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=C([CH2])C(=[CH])[CH2](16900)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=[CH](9646)', '[CH]=C([CH2])[CH2](16419)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[C]=[C]CC([CH])=C(17889)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C([CH2])[CH][C]=[CH](17890)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH][C]1CC(=[CH])C1(17891)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH][C]1C[C]=CC1(17892)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]C1([CH2])CC1=[CH](17893)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(117.479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 115.8 to 117.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]C1([CH2])C=[C]C1(17875)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(80.3914,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 77.3 to 80.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([CH2])=CC#C(17865)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C3H3(5450)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.0130473,'m^3/(mol*s)'), n=2.32138, Ea=(27.0133,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C][CH2](16918)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH]C([CH2])=C[C]=[CH](17867)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]C(=[CH])C[C]=[CH](17894)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.63133e+13,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]C([CH2])=CC=[CH](17869)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18474e+10,'s^-1'), n=0.962222, Ea=(117.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=C([CH2])C=[C][CH2](16906)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(=[CH])CC=[CH](17895)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(=[CH])C[C]=C(17482)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.912e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC[C]=[CH](16902)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.06e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(T)(20)', '[CH]=[C]C[C]=[CH](17536)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[C]=C([CH2])C[C]=[CH](17896)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=[C]CC1=CC1(17897)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=C1CC(=[CH])C1(17898)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]C1=CC(=C)C1(17860)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=C=CC(=[CH])C(17899)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=C=CC([CH2])=C(17497)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(T)(20)', '[CH]=[C]CC#C(17535)'],
    products = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]C(C)=C[C]=[CH](17871)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]CC(=[CH])[CH2](16901)'],
    products = ['[CH]=[C][CH]C([CH2])=C(17498)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4123',
    isomers = [
        '[CH]=[C]CC(=[CH])[CH2](16901)',
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
    label = 'PDepNetwork #4123',
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


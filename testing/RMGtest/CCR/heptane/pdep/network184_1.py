species(
    label = '[CH]C(=C)C[CH2](287)',
    structure = SMILES('[CH]C(=C)C[CH2]'),
    E0 = (520.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,524.801,524.987,525.161,525.162],'cm^-1')),
        HinderedRotor(inertia=(0.26622,'amu*angstrom^2'), symmetry=1, barrier=(52.167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266605,'amu*angstrom^2'), symmetry=1, barrier=(52.1713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266436,'amu*angstrom^2'), symmetry=1, barrier=(52.1704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80742,0.0413074,-1.07219e-05,-9.40701e-09,5.01172e-12,62668.9,21.6031], Tmin=(100,'K'), Tmax=(1091.36,'K')), NASAPolynomial(coeffs=[8.80615,0.0287692,-1.15122e-05,2.08538e-09,-1.42871e-13,60360.3,-16.3514], Tmin=(1091.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(17)',
    structure = SMILES('C=C'),
    E0 = (42.1493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97974,-0.0075756,5.52973e-05,-6.36222e-08,2.31767e-11,5077.46,4.04623], Tmin=(100,'K'), Tmax=(940.446,'K')), NASAPolynomial(coeffs=[5.20299,0.00782443,-2.12684e-06,3.79691e-10,-2.94671e-14,3936.28,-6.6241], Tmin=(940.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.1493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=C=C(264)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,750,770,3400,2100,316.125],'cm^-1')),
        HinderedRotor(inertia=(0.00168831,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62047e-06,-1.82971e-08,7.68641e-12,39535.3,8.27852], Tmin=(100,'K'), Tmax=(960.557,'K')), NASAPolynomial(coeffs=[6.38512,0.00814484,-2.78733e-06,4.95345e-10,-3.50146e-14,38483.6,-8.79389], Tmin=(960.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
)

species(
    label = 'CH2(T)(15)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.98,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C([CH2])=C(754)',
    structure = SMILES('[CH]C([CH2])=C'),
    E0 = (489.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,336.248,337.699,348.653],'cm^-1')),
        HinderedRotor(inertia=(0.621644,'amu*angstrom^2'), symmetry=1, barrier=(50.7924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614952,'amu*angstrom^2'), symmetry=1, barrier=(50.7076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48525,0.02547,1.30544e-05,-3.21043e-08,1.33808e-11,58905,14.749], Tmin=(100,'K'), Tmax=(981.432,'K')), NASAPolynomial(coeffs=[8.34337,0.0203035,-7.64405e-06,1.37961e-09,-9.63621e-14,56854.1,-17.9938], Tmin=(981.432,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'H(9)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]CC([CH])=C(755)',
    structure = SMILES('[CH]CC([CH])=C'),
    E0 = (763.323,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71731,0.0436494,-1.99321e-05,-1.67054e-09,2.87964e-12,91894.2,21.1432], Tmin=(100,'K'), Tmax=(1086.15,'K')), NASAPolynomial(coeffs=[9.80857,0.0251484,-9.98297e-06,1.79859e-09,-1.22908e-13,89470.2,-21.6287], Tmin=(1086.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[CH][CH2](756)',
    structure = SMILES('[CH]C([CH2])[CH][CH2]'),
    E0 = (879.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,241.787,242.702,1565.85,1566.03],'cm^-1')),
        HinderedRotor(inertia=(6.8795e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00561863,'amu*angstrom^2'), symmetry=1, barrier=(9.78223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0386542,'amu*angstrom^2'), symmetry=1, barrier=(67.2912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0028754,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9938,0.0435135,-3.26328e-05,1.41523e-08,-2.65241e-12,105829,25.8975], Tmin=(100,'K'), Tmax=(1221.13,'K')), NASAPolynomial(coeffs=[7.62383,0.0250715,-9.97927e-06,1.78472e-09,-1.20429e-13,104454,-2.38885], Tmin=(1221.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(879.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1CCC1(757)',
    structure = SMILES('[CH][C]1CCC1'),
    E0 = (614.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62891,0.0208918,3.07855e-05,-4.67252e-08,1.70605e-11,73987.4,18.1983], Tmin=(100,'K'), Tmax=(1024.83,'K')), NASAPolynomial(coeffs=[7.11543,0.0262744,-1.06014e-05,1.99525e-09,-1.41874e-13,71865.6,-9.42218], Tmin=(1024.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])CC1(758)',
    structure = SMILES('[CH]C1([CH2])CC1'),
    E0 = (626.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14118,0.0245945,4.37015e-05,-7.70622e-08,3.18216e-11,75447.1,18.5114], Tmin=(100,'K'), Tmax=(958.485,'K')), NASAPolynomial(coeffs=[14.2488,0.0147801,-4.65529e-06,8.89163e-10,-6.95211e-14,71255.9,-49.1436], Tmin=(958.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Neopentyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C=C(626)',
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
    label = '[CH]=[C][CH2](751)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.26,228.889,229.072],'cm^-1')),
        HinderedRotor(inertia=(1.35217,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27542,0.0127953,9.49555e-06,-1.56031e-08,5.42959e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.3,'K')), NASAPolynomial(coeffs=[4.1896,0.0168436,-6.77768e-06,1.2222e-09,-8.33566e-14,73336.4,4.8934], Tmin=(1063.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][CH2](18)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (313.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1608.02,1609.59,1610.81,1612.2,1612.86],'cm^-1')),
        HinderedRotor(inertia=(0.00575811,'amu*angstrom^2'), symmetry=1, barrier=(10.593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86546,-0.00370153,4.71859e-05,-6.15022e-08,2.54427e-11,37752.3,8.63154], Tmin=(100,'K'), Tmax=(845.133,'K')), NASAPolynomial(coeffs=[5.8996,0.00441368,1.29131e-06,-4.57987e-10,3.67865e-14,36774.8,-4.58851], Tmin=(845.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ)"""),
)

species(
    label = '[CH]C(=C)[CH][CH2](759)',
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
    label = '[CH]C(=[CH])C[CH2](760)',
    structure = SMILES('[CH]C(=[CH])C[CH2]'),
    E0 = (767.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,499.609,499.614,499.619],'cm^-1')),
        HinderedRotor(inertia=(0.294249,'amu*angstrom^2'), symmetry=1, barrier=(52.1221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29425,'amu*angstrom^2'), symmetry=1, barrier=(52.1219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294248,'amu*angstrom^2'), symmetry=1, barrier=(52.1221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65878,0.0455874,-2.94386e-05,9.81929e-09,-1.35238e-12,92392,21.8548], Tmin=(100,'K'), Tmax=(1655.31,'K')), NASAPolynomial(coeffs=[11.0697,0.0228462,-8.83111e-06,1.51975e-09,-9.89076e-14,89276.4,-28.2906], Tmin=(1655.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)[CH]C(761)',
    structure = SMILES('[CH]C([CH2])=CC'),
    E0 = (453.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77232,0.0417854,-9.86968e-06,-1.09886e-08,5.66972e-12,54596.9,19.1164], Tmin=(100,'K'), Tmax=(1079.7,'K')), NASAPolynomial(coeffs=[8.99265,0.029042,-1.16235e-05,2.10866e-09,-1.44757e-13,52221.3,-20.0517], Tmin=(1079.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])CC(762)',
    structure = SMILES('[CH]C(=[CH])CC'),
    E0 = (562.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,366.525,366.526,366.527],'cm^-1')),
        HinderedRotor(inertia=(0.538496,'amu*angstrom^2'), symmetry=1, barrier=(51.336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538505,'amu*angstrom^2'), symmetry=1, barrier=(51.336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538507,'amu*angstrom^2'), symmetry=1, barrier=(51.336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68232,0.0442599,-1.74826e-05,-3.58511e-09,3.25014e-12,67706.6,19.9758], Tmin=(100,'K'), Tmax=(1106.97,'K')), NASAPolynomial(coeffs=[9.35354,0.0280645,-1.1153e-05,2.0075e-09,-1.36834e-13,65302.1,-21.0023], Tmin=(1106.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC[CH2](280)',
    structure = SMILES('[CH]=[C]CC[CH2]'),
    E0 = (649.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.38393,'amu*angstrom^2'), symmetry=1, barrier=(31.8193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000321969,'amu*angstrom^2'), symmetry=1, barrier=(3.65565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714807,'amu*angstrom^2'), symmetry=1, barrier=(16.4348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3052.11,'J/mol'), sigma=(5.53315,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.73 K, Pc=40.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07344,0.0423233,-2.94338e-05,1.10648e-08,-1.78046e-12,78231.9,21.8978], Tmin=(100,'K'), Tmax=(1387.7,'K')), NASAPolynomial(coeffs=[8.19337,0.0246827,-1.03656e-05,1.90418e-09,-1.30139e-13,76533.4,-9.63256], Tmin=(1387.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C[CH2](308)',
    structure = SMILES('[CH]=[C]C[CH2]'),
    E0 = (673.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.277912,'amu*angstrom^2'), symmetry=1, barrier=(6.38974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.020326,'amu*angstrom^2'), symmetry=1, barrier=(9.87347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6422,0.0283978,-1.83267e-05,6.21845e-09,-8.8654e-13,81072.8,17.6148], Tmin=(100,'K'), Tmax=(1566.75,'K')), NASAPolynomial(coeffs=[7.52448,0.015933,-6.39296e-06,1.14053e-09,-7.62759e-14,79542.9,-8.13161], Tmin=(1566.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[C]=C([CH2])C[CH2](763)',
    structure = SMILES('[C]C(=C)C[CH2]'),
    E0 = (819.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,797.162,798.794],'cm^-1')),
        HinderedRotor(inertia=(0.600105,'amu*angstrom^2'), symmetry=1, barrier=(13.7976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600572,'amu*angstrom^2'), symmetry=1, barrier=(13.8083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75745,0.0434769,-3.12037e-05,9.53745e-09,-5.40046e-13,98606.3,19.1262], Tmin=(100,'K'), Tmax=(1115.15,'K')), NASAPolynomial(coeffs=[11.1626,0.0175433,-6.81548e-06,1.23182e-09,-8.46579e-14,96023.6,-29.4487], Tmin=(1115.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(819.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CJ3)"""),
)

species(
    label = '[CH2]CC1=CC1(764)',
    structure = SMILES('[CH2]CC1=CC1'),
    E0 = (418.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16556,0.0296296,1.69712e-05,-4.23206e-08,1.80755e-11,50373.1,18.7945], Tmin=(100,'K'), Tmax=(975.312,'K')), NASAPolynomial(coeffs=[11.0442,0.0191736,-6.86947e-06,1.26365e-09,-9.1096e-14,47406.7,-30.1469], Tmin=(975.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C1CCC1(748)',
    structure = SMILES('[CH]=C1CCC1'),
    E0 = (363.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62009,0.0185523,4.26385e-05,-6.24352e-08,2.33412e-11,43732.2,15.6225], Tmin=(100,'K'), Tmax=(998.379,'K')), NASAPolynomial(coeffs=[8.6752,0.0237138,-9.31995e-06,1.77694e-09,-1.29063e-13,41056.8,-20.9235], Tmin=(998.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=CCC1(765)',
    structure = SMILES('C=C1[CH]CC1'),
    E0 = (257.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90085,0.0097037,6.85862e-05,-8.87071e-08,3.23833e-11,30977.6,13.798], Tmin=(100,'K'), Tmax=(987.001,'K')), NASAPolynomial(coeffs=[8.44602,0.0241142,-9.36782e-06,1.80785e-09,-1.33445e-13,28086.5,-21.9826], Tmin=(987.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=C(C)C=C(561)',
    structure = SMILES('[CH]=C(C)C=C'),
    E0 = (303.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.896726,'amu*angstrom^2'), symmetry=1, barrier=(20.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895838,'amu*angstrom^2'), symmetry=1, barrier=(20.5971,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76567,0.0404503,-9.85362e-06,-1.81195e-08,1.05257e-11,36582.1,17.0183], Tmin=(100,'K'), Tmax=(967.447,'K')), NASAPolynomial(coeffs=[12.2564,0.0175176,-5.99196e-06,1.06038e-09,-7.45689e-14,33595.6,-38.1904], Tmin=(967.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=C)C=C(292)',
    structure = SMILES('[CH2]C(=C)C=C'),
    E0 = (207.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.21475,'amu*angstrom^2'), symmetry=1, barrier=(27.9294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21659,'amu*angstrom^2'), symmetry=1, barrier=(27.9719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85253,0.0342665,1.60034e-05,-4.98275e-08,2.29338e-11,25085.3,15.8911], Tmin=(100,'K'), Tmax=(944.855,'K')), NASAPolynomial(coeffs=[14.18,0.0145156,-4.13607e-06,7.15936e-10,-5.30912e-14,21307.9,-50.5445], Tmin=(944.855,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC[CH2](306)',
    structure = SMILES('C#CC[CH2]'),
    E0 = (354.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.0191618,'amu*angstrom^2'), symmetry=1, barrier=(13.4063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.85815,'amu*angstrom^2'), symmetry=1, barrier=(65.7146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72207,0.023801,-9.48739e-07,-1.56804e-08,8.25241e-12,42725.9,14.2949], Tmin=(100,'K'), Tmax=(909.863,'K')), NASAPolynomial(coeffs=[7.43485,0.0145133,-4.48185e-06,7.16223e-10,-4.67695e-14,41395.2,-10.5966], Tmin=(909.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C(C)[CH][CH2](562)',
    structure = SMILES('[CH]C(C)=C[CH2]'),
    E0 = (453.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,361.815,361.816,361.816],'cm^-1')),
        HinderedRotor(inertia=(0.549955,'amu*angstrom^2'), symmetry=1, barrier=(51.0891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549956,'amu*angstrom^2'), symmetry=1, barrier=(51.0891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549956,'amu*angstrom^2'), symmetry=1, barrier=(51.0891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77232,0.0417854,-9.86968e-06,-1.09886e-08,5.66972e-12,54596.9,19.1164], Tmin=(100,'K'), Tmax=(1079.7,'K')), NASAPolynomial(coeffs=[8.99265,0.029042,-1.16235e-05,2.10866e-09,-1.44757e-13,52221.3,-20.0517], Tmin=(1079.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])=C(293)',
    structure = SMILES('[CH2]C=C([CH2])[CH2]'),
    E0 = (385.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.110043,'amu*angstrom^2'), symmetry=1, barrier=(31.0187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110116,'amu*angstrom^2'), symmetry=1, barrier=(31.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242242,'amu*angstrom^2'), symmetry=1, barrier=(68.3077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82336,0.036309,8.59652e-06,-3.90391e-08,1.80728e-11,46458.6,18.3828], Tmin=(100,'K'), Tmax=(965.965,'K')), NASAPolynomial(coeffs=[13.046,0.0178172,-6.13813e-06,1.11719e-09,-8.09237e-14,42985,-42.1286], Tmin=(965.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
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
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (520.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (927.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (975.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (902.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (646.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (626.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (648.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (678.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (928.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (816.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (979.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (634.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (606.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (895.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1111.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1030.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (528.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (528.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (528.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (583.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (583.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (752.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (666.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (722.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (712.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['C=C(17)', '[CH]=C=C(264)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(15)', '[CH]C([CH2])=C(754)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.45208e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(9)', '[CH]CC([CH])=C(755)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C([CH2])[CH][CH2](756)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH][C]1CCC1(757)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.37682e+08,'s^-1'), n=1.00167, Ea=(126.008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH]C1([CH2])CC1(758)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(106.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 103.9 to 106.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(9)', '[CH]=C([CH2])C=C(626)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+08,'cm^3/(mol*s)'), n=1.64, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2580 used for Cds-CdH_Cds-HH;HJ
Exact match found for rate rule [Cds-CdH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C(17)', '[CH]=[C][CH2](751)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00310793,'m^3/(mol*s)'), n=2.49201, Ea=(21.4614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH2](18)', '[CH]=[C][CH2](751)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(9)', '[CH]C(=C)[CH][CH2](759)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(9)', '[CH]C(=[CH])C[CH2](760)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.63133e+13,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH]C(=C)[CH]C(761)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(=[CH])CC(762)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(222600,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2(T)(15)', '[CH]=[C]C[CH2](308)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(9)', '[C]=C([CH2])C[CH2](763)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH2]CC1=CC1(764)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH]=C1CCC1(748)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH2]C1=CCC1(765)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH]=C(C)C=C(561)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH2]C(=C)C=C(292)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(T)(15)', 'C#CC[CH2](306)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH2](18)', '[CH]=C=C(264)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00625461,'m^3/(mol*s)'), n=2.54618, Ea=(24.6367,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cdd_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH]=C(C)[CH][CH2](562)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]C(=C)C[CH2](287)'],
    products = ['[CH2][CH]C([CH2])=C(293)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #184',
    isomers = [
        '[CH]C(=C)C[CH2](287)',
    ],
    reactants = [
        ('C=C(17)', '[CH]=C=C(264)'),
    ],
    bathGas = {
        'Ar': 0.25,
        'He': 0.25,
        'Ne': 0.25,
        'N2': 0.25,
    },
)

pressureDependence(
    label = 'PDepNetwork #184',
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


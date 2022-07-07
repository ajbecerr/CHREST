species(
    label = '[CH2]C(=C)C([CH2])([O])[CH][O](25001)',
    structure = SMILES('[CH2]C(=C)C([CH2])([O])[CH][O]'),
    E0 = (567.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0485096,0.0940914,-0.000137715,1.14908e-07,-3.8559e-11,68353.9,33.6046], Tmin=(100,'K'), Tmax=(815.343,'K')), NASAPolynomial(coeffs=[9.6972,0.0378146,-1.77322e-05,3.3537e-09,-2.30446e-13,67077.8,-9.15242], Tmin=(815.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4007.73,'J/mol'), sigma=(6.4029,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=626.00 K, Pc=34.64 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'allene(458)',
    structure = SMILES('C=C=C'),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.0070462,2.78306e-05,-3.99445e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.705,'K')), NASAPolynomial(coeffs=[6.79956,0.00959979,-3.02068e-06,5.37827e-10,-3.92606e-14,19772.3,-12.7582], Tmin=(949.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(=C)C([CH2])([O])C=O(25013)',
    structure = SMILES('[CH2]C(=C)C([CH2])([O])C=O'),
    E0 = (265.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183147,0.100269,-0.000155673,1.31978e-07,-4.39765e-11,32100.1,29.8889], Tmin=(100,'K'), Tmax=(846.83,'K')), NASAPolynomial(coeffs=[10.5635,0.0356189,-1.65573e-05,3.09261e-09,-2.09825e-13,30778,-17.2303], Tmin=(846.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=CC(O)(C=O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C[C]([O])[CH][O](25974)',
    structure = SMILES('[CH2]C(=C)C[C]([O])[CH][O]'),
    E0 = (552.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,360,370,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4370.9,'J/mol'), sigma=(7.31103,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.73 K, Pc=25.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107255,0.0996473,-0.000156385,1.36358e-07,-4.68006e-11,66629.7,32.1127], Tmin=(100,'K'), Tmax=(838.265,'K')), NASAPolynomial(coeffs=[9.16542,0.0390369,-1.86476e-05,3.52891e-09,-2.41389e-13,65650.1,-7.55752], Tmin=(838.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(Allyl_P) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]([O])C([O])C([CH2])=C(25956)',
    structure = SMILES('[CH2][C]([O])C([O])C([CH2])=C'),
    E0 = (575.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,180,876.961,877.024],'cm^-1')),
        HinderedRotor(inertia=(0.033989,'amu*angstrom^2'), symmetry=1, barrier=(18.5514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806874,'amu*angstrom^2'), symmetry=1, barrier=(18.5516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00809045,'amu*angstrom^2'), symmetry=1, barrier=(4.41567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.806879,'amu*angstrom^2'), symmetry=1, barrier=(18.5517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0533025,0.0883307,-9.90476e-05,5.68824e-08,-1.29371e-11,69399.9,32.5192], Tmin=(100,'K'), Tmax=(1071.41,'K')), NASAPolynomial(coeffs=[16.804,0.0253948,-1.09344e-05,2.05461e-09,-1.43512e-13,65787.8,-49.9698], Tmin=(1071.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])([CH][O])C[C]=C(24956)',
    structure = SMILES('[CH2]C([O])([CH][O])C[C]=C'),
    E0 = (670.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180,945.309,1291.46,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161898,'amu*angstrom^2'), symmetry=1, barrier=(3.72236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161898,'amu*angstrom^2'), symmetry=1, barrier=(3.72236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161898,'amu*angstrom^2'), symmetry=1, barrier=(3.72236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161898,'amu*angstrom^2'), symmetry=1, barrier=(3.72236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4372.6,'J/mol'), sigma=(7.32671,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.99 K, Pc=25.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.379689,0.109135,-0.000186942,1.6975e-07,-5.90408e-11,80777.5,33.8989], Tmin=(100,'K'), Tmax=(861.112,'K')), NASAPolynomial(coeffs=[8.68495,0.0400147,-1.94837e-05,3.67521e-09,-2.49097e-13,80217.9,-2.66198], Tmin=(861.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S)"""),
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
    label = '[CH2][C]([CH][O])C([CH2])=C(27802)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])[CH][O]'),
    E0 = (536.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,570.948,571.378],'cm^-1')),
        HinderedRotor(inertia=(0.00051542,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000517073,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211605,'amu*angstrom^2'), symmetry=1, barrier=(48.9508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21107,'amu*angstrom^2'), symmetry=1, barrier=(48.9323,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864074,0.0586113,-3.00475e-05,-3.47792e-09,5.45461e-12,64629.8,27.7985], Tmin=(100,'K'), Tmax=(1049.19,'K')), NASAPolynomial(coeffs=[14.9996,0.0246455,-9.97446e-06,1.86822e-09,-1.32275e-13,60567,-46.3026], Tmin=(1049.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(CCOJ) + radical(Allyl_P) + radical(C=CCJO) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2][C]([CH2])C([O])=C[O](12296)',
    structure = SMILES('[CH2]C([CH2])=C([O])[CH][O]'),
    E0 = (347.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,422.949,422.954,422.973],'cm^-1')),
        HinderedRotor(inertia=(0.000942335,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000942187,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000942311,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07031,0.0559424,-3.87043e-05,5.9337e-09,2.81786e-12,41865.7,27.5769], Tmin=(100,'K'), Tmax=(1010.07,'K')), NASAPolynomial(coeffs=[15.0054,0.0175457,-6.61405e-06,1.20824e-09,-8.52133e-14,38194.2,-44.0307], Tmin=(1010.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH][O](1548)',
    structure = SMILES('[CH][O]'),
    E0 = (431.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([806.798,806.798,2696.84],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86882,-0.000695042,1.53508e-05,-2.00576e-08,7.65929e-12,51865.1,6.76618], Tmin=(100,'K'), Tmax=(973.285,'K')), NASAPolynomial(coeffs=[6.03251,-0.000301449,4.32946e-07,-3.67045e-11,-1.25373e-15,51004.1,-5.87326], Tmin=(973.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[CH2][C]([CH2])C(=C)[O](2601)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])[O]'),
    E0 = (316.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.233452,'amu*angstrom^2'), symmetry=1, barrier=(19.9043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0183169,'amu*angstrom^2'), symmetry=1, barrier=(20.0446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0723988,'amu*angstrom^2'), symmetry=1, barrier=(78.9784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.504591,0.0616067,-5.89295e-05,2.8048e-08,-5.05278e-12,38169.5,25.0932], Tmin=(100,'K'), Tmax=(1549.14,'K')), NASAPolynomial(coeffs=[16.8162,0.0112712,-2.23363e-06,2.24815e-10,-1.00569e-14,34101.8,-57.5577], Tmin=(1549.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])([C]=C)[CH][O](27024)',
    structure = SMILES('[CH2]C([O])([C]=C)[CH][O]'),
    E0 = (688.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180,180,180,232.238,945.888,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.147625,'amu*angstrom^2'), symmetry=1, barrier=(3.39418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147625,'amu*angstrom^2'), symmetry=1, barrier=(3.39418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147625,'amu*angstrom^2'), symmetry=1, barrier=(3.39418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778114,0.0821809,-0.000144161,1.354e-07,-4.82839e-11,82902.3,31.1128], Tmin=(100,'K'), Tmax=(864.004,'K')), NASAPolynomial(coeffs=[5.71281,0.0340421,-1.66753e-05,3.14949e-09,-2.13447e-13,82993.7,13.4904], Tmin=(864.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(Cds_S)"""),
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
    label = '[CH]C([O])([CH][O])C([CH2])=C(27803)',
    structure = SMILES('[CH]C([O])([CH][O])C([CH2])=C'),
    E0 = (801.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.136268,0.0909467,-0.0001307,1.04889e-07,-3.39719e-11,96588.5,33.4204], Tmin=(100,'K'), Tmax=(793.511,'K')), NASAPolynomial(coeffs=[10.9005,0.0333787,-1.56269e-05,2.95976e-09,-2.03935e-13,94984.3,-15.3654], Tmin=(793.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(801.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])([O])[C][O](27804)',
    structure = SMILES('[CH2]C(=C)C([CH2])([O])[C][O]'),
    E0 = (847.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.278185,0.0888873,-0.000121168,8.84781e-08,-2.49826e-11,102105,31.9851], Tmin=(100,'K'), Tmax=(674.865,'K')), NASAPolynomial(coeffs=[11.1138,0.0331067,-1.5954e-05,3.08114e-09,-2.15515e-13,100451,-17.4542], Tmin=(674.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(Allyl_P) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C(=C)C([CH2])([O])[CH][O](25158)',
    structure = SMILES('[CH]C(=C)C([CH2])([O])[CH][O]'),
    E0 = (786.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,180,1096.77,1096.77,1096.77,1096.77,1096.77,1096.77,1096.77,1096.77,1096.77,1096.77,2269.21],'cm^-1')),
        HinderedRotor(inertia=(0.0953526,'amu*angstrom^2'), symmetry=1, barrier=(2.19234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953526,'amu*angstrom^2'), symmetry=1, barrier=(2.19234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953526,'amu*angstrom^2'), symmetry=1, barrier=(2.19234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953526,'amu*angstrom^2'), symmetry=1, barrier=(2.19234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0326055,0.0967193,-0.000146702,1.28875e-07,-4.4862e-11,94714.2,35.1818], Tmin=(100,'K'), Tmax=(842.177,'K')), NASAPolynomial(coeffs=[7.38616,0.0438856,-2.07058e-05,3.88878e-09,-2.64978e-13,94110.6,4.73853], Tmin=(842.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(786.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)C1([CH][O])CO1(27805)',
    structure = SMILES('[CH2]C(=C)C1([CH][O])CO1'),
    E0 = (310.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.366994,0.0780372,-7.13358e-05,2.32709e-08,3.27804e-12,37421.6,27.431], Tmin=(100,'K'), Tmax=(784.655,'K')), NASAPolynomial(coeffs=[15.0094,0.0236574,-6.11767e-06,7.72347e-10,-4.00568e-14,34499.9,-43.6341], Tmin=(784.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1([CH2])OC1[O](24187)',
    structure = SMILES('[CH2]C(=C)C1([CH2])OC1[O]'),
    E0 = (299.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4117.5,'J/mol'), sigma=(6.93639,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=643.14 K, Pc=27.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26811,0.0781585,-6.92136e-05,2.06602e-08,3.57149e-12,36201.2,27.7447], Tmin=(100,'K'), Tmax=(821.339,'K')), NASAPolynomial(coeffs=[16.2773,0.0216813,-5.31496e-06,6.49292e-10,-3.34487e-14,32846.6,-50.7517], Tmin=(821.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1([O])CC1[O](27806)',
    structure = SMILES('[CH2]C(=C)C1([O])CC1[O]'),
    E0 = (323.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687434,0.0568864,-4.88113e-06,-4.10816e-08,2.15181e-11,39088,28.3474], Tmin=(100,'K'), Tmax=(965.184,'K')), NASAPolynomial(coeffs=[18.9206,0.0196422,-6.55171e-06,1.20566e-09,-8.93232e-14,33783.5,-68.217], Tmin=(965.184,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CC(C)2OJ) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([CH][O])OCC1=C(27715)',
    structure = SMILES('[CH2]C1([CH][O])OCC1=C'),
    E0 = (375.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.692615,0.0732903,-6.97818e-05,3.62496e-08,-7.76639e-12,45233.9,27.9762], Tmin=(100,'K'), Tmax=(1109.6,'K')), NASAPolynomial(coeffs=[12.1944,0.0318278,-1.37314e-05,2.57363e-09,-1.79009e-13,42681.4,-28.7092], Tmin=(1109.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = 'C=C1CCC1([O])[CH][O](27629)',
    structure = SMILES('C=C1CCC1([O])[CH][O]'),
    E0 = (369.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51414,0.0600141,-4.28221e-05,1.67478e-08,-2.93445e-12,44545.3,26.6176], Tmin=(100,'K'), Tmax=(1227.35,'K')), NASAPolynomial(coeffs=[7.60283,0.0401709,-1.85709e-05,3.57528e-09,-2.51331e-13,43050.7,-4.00409], Tmin=(1227.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])C(=C)CC1[O](27781)',
    structure = SMILES('[CH2]C1([O])C(=C)CC1[O]'),
    E0 = (394.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719031,0.0624474,-3.78849e-05,6.25183e-09,1.4616e-12,47612.6,28.7127], Tmin=(100,'K'), Tmax=(1134.08,'K')), NASAPolynomial(coeffs=[14.7235,0.0277495,-1.1431e-05,2.12846e-09,-1.48565e-13,43491,-44.7797], Tmin=(1134.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)2OJ) + radical(CC(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([O])([CH][O])[C]1CC1(27807)',
    structure = SMILES('[CH2]C([O])([CH][O])[C]1CC1'),
    E0 = (610.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.606111,0.0776417,-7.98558e-05,4.522e-08,-1.06034e-11,73596.2,28.8937], Tmin=(100,'K'), Tmax=(1016.51,'K')), NASAPolynomial(coeffs=[11.706,0.0339636,-1.54032e-05,2.94971e-09,-2.07491e-13,71339.6,-24.8385], Tmin=(1016.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CJC(C)2O) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]1COC1([CH2])[CH][O](27808)',
    structure = SMILES('[CH2][C]1COC1([CH2])[CH][O]'),
    E0 = (607.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634845,0.0776763,-7.22382e-05,1.52548e-08,1.38605e-11,73125,28.3311], Tmin=(100,'K'), Tmax=(661.816,'K')), NASAPolynomial(coeffs=[12.8007,0.026559,-7.179e-06,8.88913e-10,-4.22685e-14,71023.8,-29.0485], Tmin=(661.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CJC(C)OC) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1CCC1([O])[CH][O](27809)',
    structure = SMILES('[CH2][C]1CCC1([O])[CH][O]'),
    E0 = (599.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732658,0.0655561,-4.98675e-05,1.99782e-08,-3.27815e-12,72233.7,30.3905], Tmin=(100,'K'), Tmax=(1428.19,'K')), NASAPolynomial(coeffs=[13.8053,0.0289427,-1.14131e-05,2.02799e-09,-1.36007e-13,68499.7,-37.3369], Tmin=(1428.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1CC([O])C1([CH2])[O](27810)',
    structure = SMILES('[CH2][C]1CC([O])C1([CH2])[O]'),
    E0 = (623.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5266,0.0636498,-2.58022e-05,-1.79808e-08,1.31775e-11,75101.7,29.3228], Tmin=(100,'K'), Tmax=(965.609,'K')), NASAPolynomial(coeffs=[17.756,0.0220697,-7.48997e-06,1.328e-09,-9.41168e-14,70385.4,-60.3878], Tmin=(965.609,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(CJC(C)2O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1CO[CH]C1([CH2])[O](27811)',
    structure = SMILES('[CH2][C]1CO[CH]C1([CH2])[O]'),
    E0 = (533.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.648517,0.0795032,-7.6115e-05,3.79516e-08,-7.06327e-12,64343.3,29.8586], Tmin=(100,'K'), Tmax=(1581.73,'K')), NASAPolynomial(coeffs=[16.9542,0.0180805,-1.83188e-06,-1.15305e-10,2.15382e-14,60889.8,-56.4496], Tmin=(1581.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(533.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)2OJ) + radical(CCJ(C)CO) + radical(CCsJOCs) + radical(CJC(C)2O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([CH2])OC1([CH2])[CH][O](27735)',
    structure = SMILES('[CH2]C1([CH2])OC1([CH2])[CH][O]'),
    E0 = (654.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53874,0.133411,-0.000232787,1.93389e-07,-5.90268e-11,78886,29.5416], Tmin=(100,'K'), Tmax=(983.036,'K')), NASAPolynomial(coeffs=[16.1227,0.0234681,-6.92321e-06,8.10098e-10,-3.03948e-14,77253.5,-46.0046], Tmin=(983.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([CH2])CC1([O])[CH][O](27649)',
    structure = SMILES('[CH2]C1([CH2])CC1([O])[CH][O]'),
    E0 = (644.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493579,0.0792079,-8.26549e-05,4.61979e-08,-1.05427e-11,77611.6,29.9963], Tmin=(100,'K'), Tmax=(1050.06,'K')), NASAPolynomial(coeffs=[13.0793,0.0312663,-1.41724e-05,2.72051e-09,-1.91815e-13,74968.4,-31.3375], Tmin=(1050.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(Neopentyl) + radical(Neopentyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([CH2])C([O])C1([CH2])[O](27812)',
    structure = SMILES('[CH2]C1([CH2])C([O])C1([CH2])[O]'),
    E0 = (667.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.392464,0.0852375,-8.57922e-05,4.24544e-08,-8.1376e-12,80509,31.3714], Tmin=(100,'K'), Tmax=(1282.59,'K')), NASAPolynomial(coeffs=[21.4316,0.0171753,-6.1929e-06,1.0802e-09,-7.30271e-14,74910.7,-79.3485], Tmin=(1282.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(Neopentyl) + radical(Neopentyl) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])[CH]OC1([CH2])[CH2](27685)',
    structure = SMILES('[CH2]C1([O])[CH]OC1([CH2])[CH2]'),
    E0 = (653.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.41626,0.117156,-0.000170059,1.19071e-07,-3.09495e-11,78814.9,29.3379], Tmin=(100,'K'), Tmax=(1097.74,'K')), NASAPolynomial(coeffs=[21.4836,0.0135792,-1.01686e-06,-2.96863e-10,4.03897e-14,75000.3,-77.7509], Tmin=(1097.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)2OJ) + radical(CCsJOCs) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C(=C)C([O])=C[O](12294)',
    structure = SMILES('[CH2]C(=C)C([O])=C[O]'),
    E0 = (64.4117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30302,'amu*angstrom^2'), symmetry=1, barrier=(29.9591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29734,'amu*angstrom^2'), symmetry=1, barrier=(29.8284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.70548,0.0763064,-8.01811e-05,3.90501e-08,-6.92726e-12,7939.82,27.8012], Tmin=(100,'K'), Tmax=(1656,'K')), NASAPolynomial(coeffs=[21.9763,0.00356412,1.97243e-06,-5.70396e-10,4.25168e-14,2889.6,-85.6326], Tmin=(1656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C(=C)[O](2600)',
    structure = SMILES('[CH2]C(=C)C(=C)[O]'),
    E0 = (131.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01442,'amu*angstrom^2'), symmetry=1, barrier=(23.3234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01583,'amu*angstrom^2'), symmetry=1, barrier=(23.3559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.45,'J/mol'), sigma=(6.08582,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.10 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15813,0.0532368,-3.04307e-05,-7.62019e-09,9.61083e-12,15955.7,19.1278], Tmin=(100,'K'), Tmax=(927.538,'K')), NASAPolynomial(coeffs=[16.1032,0.012595,-3.20818e-06,4.87529e-10,-3.33958e-14,12159.2,-57.3701], Tmin=(927.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C(=C)[CH][O](23578)',
    structure = SMILES('[CH2]C(=C)C([CH2])=C[O]'),
    E0 = (254.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43333,'amu*angstrom^2'), symmetry=1, barrier=(32.9551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4326,'amu*angstrom^2'), symmetry=1, barrier=(32.9383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43095,'amu*angstrom^2'), symmetry=1, barrier=(32.9003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40155,0.0603326,-6.20328e-06,-5.23668e-08,2.94114e-11,30745.2,23.9436], Tmin=(100,'K'), Tmax=(925.932,'K')), NASAPolynomial(coeffs=[23.7407,0.00806211,-1.82938e-07,-6.85501e-11,2.11033e-16,24341.7,-98.098], Tmin=(925.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=C(459)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.110742,'amu*angstrom^2'), symmetry=1, barrier=(23.6192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28575,0.0103679,1.6796e-05,-2.70041e-08,1.05227e-11,47575.3,10.4038], Tmin=(100,'K'), Tmax=(980.899,'K')), NASAPolynomial(coeffs=[6.52804,0.0104286,-3.60844e-06,6.68464e-10,-4.85509e-14,46300.2,-8.43297], Tmin=(980.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])[CH][O](26240)',
    structure = SMILES('[CH2][C]([O])[CH][O]'),
    E0 = (571.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,180,2035.01],'cm^-1')),
        HinderedRotor(inertia=(0.0866393,'amu*angstrom^2'), symmetry=1, barrier=(1.99201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086712,'amu*angstrom^2'), symmetry=1, barrier=(1.99368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60835,0.065059,-0.00013406,1.3209e-07,-4.73466e-11,68841.1,21.2809], Tmin=(100,'K'), Tmax=(890.641,'K')), NASAPolynomial(coeffs=[4.29645,0.0226555,-1.15626e-05,2.16119e-09,-1.43083e-13,69565.2,15.3771], Tmin=(890.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=C(C)C([CH2])([O])[CH][O](27813)',
    structure = SMILES('[CH]=C(C)C([CH2])([O])[CH][O]'),
    E0 = (662.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,309.51,1461.24,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.141554,'amu*angstrom^2'), symmetry=1, barrier=(3.25459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141554,'amu*angstrom^2'), symmetry=1, barrier=(3.25459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141554,'amu*angstrom^2'), symmetry=1, barrier=(3.25459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141554,'amu*angstrom^2'), symmetry=1, barrier=(3.25459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.068468,0.100657,-0.000165027,1.48642e-07,-5.18747e-11,79852,34.838], Tmin=(100,'K'), Tmax=(853.59,'K')), NASAPolynomial(coeffs=[7.87449,0.0406337,-1.94779e-05,3.67127e-09,-2.49642e-13,79326.7,2.64118], Tmin=(853.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH2])C([CH2])(O)[CH][O](27814)',
    structure = SMILES('[CH]C(=C)C([CH2])(O)[CH][O]'),
    E0 = (557.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.351632,0.103017,-0.000149726,1.22714e-07,-4.00351e-11,67175.5,35.3485], Tmin=(100,'K'), Tmax=(852.145,'K')), NASAPolynomial(coeffs=[10.5735,0.0399176,-1.7855e-05,3.27425e-09,-2.20191e-13,65742.5,-13.0936], Tmin=(852.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C(C)([O])[CH][O](27815)',
    structure = SMILES('[CH]C(=C)C(C)([O])[CH][O]'),
    E0 = (572.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.209364,0.0900867,-0.000118698,9.743e-08,-3.32837e-11,69039.2,33.6745], Tmin=(100,'K'), Tmax=(799.654,'K')), NASAPolynomial(coeffs=[7.49529,0.045874,-2.11931e-05,3.99415e-09,-2.74852e-13,68122.3,1.70583], Tmin=(799.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C([CH2])([O])C[O](27816)',
    structure = SMILES('[CH]C(=C)C([CH2])([O])C[O]'),
    E0 = (606.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1094.28,1094.28,1094.28,1094.28,1094.28,1094.28,1094.28,1094.28,1094.28,1094.28,2270.36],'cm^-1')),
        HinderedRotor(inertia=(0.0942496,'amu*angstrom^2'), symmetry=1, barrier=(2.16698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0942496,'amu*angstrom^2'), symmetry=1, barrier=(2.16698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0942496,'amu*angstrom^2'), symmetry=1, barrier=(2.16698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0942496,'amu*angstrom^2'), symmetry=1, barrier=(2.16698,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.328204,0.0876507,-0.000114305,9.51034e-08,-3.31409e-11,73021.2,34.5735], Tmin=(100,'K'), Tmax=(795.301,'K')), NASAPolynomial(coeffs=[6.64157,0.0473095,-2.20205e-05,4.16819e-09,-2.87664e-13,72288.6,7.26859], Tmin=(795.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C([CH2])([O])[CH]O(27817)',
    structure = SMILES('[CH]C(=C)C([CH2])([O])[CH]O'),
    E0 = (560.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0416138,0.0955244,-0.000117553,7.62212e-08,-1.76277e-11,67573.9,34.6202], Tmin=(100,'K'), Tmax=(677.773,'K')), NASAPolynomial(coeffs=[11.6871,0.0381885,-1.69604e-05,3.14554e-09,-2.15088e-13,65711.1,-19.4159], Tmin=(677.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(AllylJ2_triplet)"""),
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
    E0 = (567.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (567.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (724.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (739.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (840.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1055.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (784.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (802.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1126.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1013.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1059.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (998.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (572.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (570.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (572.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (575.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (575.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (575.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (798.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (693.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (697.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (693.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (641.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (659.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (644.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (667.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (689.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (567.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (584.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (567.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (567.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (779.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (967.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (792.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (986.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (617.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (650.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (695.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['C=C([O])[CH][O](2850)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C(=C)C([CH2])([O])C=O(25013)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C(=C)C[C]([O])[CH][O](25974)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][C]([O])C([O])C([CH2])=C(25956)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2][C]([CH][O])C([CH2])=C(27802)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][C]([CH2])C([O])=C[O](12296)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][O](1548)', '[CH2][C]([CH2])C(=C)[O](2601)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(T)(20)', '[CH2]C([O])([C]=C)[CH][O](27024)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([O])([CH][O])C([CH2])=C(27803)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]C(=C)C([CH2])([O])[C][O](27804)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C(=C)C([CH2])([O])[CH][O](25158)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C(=C)C1([CH][O])CO1(27805)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C(=C)C1([O])CC1[O](27806)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.01367e+13,'s^-1'), n=0.0238333, Ea=(5.50893,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_2H;Ypri_rad_out] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C1([CH][O])OCC1=C(27715)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['C=C1CCC1([O])[CH][O](27629)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C1([O])C(=C)CC1[O](27781)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C([O])([CH][O])[C]1CC1(27807)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2][C]1COC1([CH2])[CH][O](27808)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.06838e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2][C]1CCC1([O])[CH][O](27809)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.43734e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2][C]1CC([O])C1([CH2])[O](27810)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.66736e+08,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2][C]1CO[CH]C1([CH2])[O](27811)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.10072e+11,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C1([CH2])OC1([CH2])[CH][O](27735)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(92.3853,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C1([CH2])CC1([O])[CH][O](27649)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(77.0629,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 75.6 to 77.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C1([CH2])C([O])C1([CH2])[O](27812)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(100.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 99.4 to 100.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C1([O])[CH]OC1([CH2])[CH2](27685)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(T)(20)', '[CH2]C(=C)C([O])=C[O](12294)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(120.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;Y_1centerbirad] for rate rule [CO_O;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.8 to 120.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH][O](1548)', '[CH2]C(=C)C(=C)[O](2600)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.75465,'m^3/(mol*s)'), n=1.70829, Ea=(22.0129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [CO_O;YJ] for rate rule [CO_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(4)', '[CH2]C(=C)C(=C)[CH][O](23578)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(69.7881,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 69.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C([O])[CH][O](2850)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(150.934,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 148.5 to 150.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]([O])[CH][O](26240)', 'allene(458)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C]([O])[CH][O](26240)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(C)C([CH2])([O])[CH][O](27813)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH]=C([CH2])C([CH2])(O)[CH][O](27814)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.84141e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_singleH] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([CH2])C(C)([O])[CH][O](27815)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C([CH2])C([CH2])([O])C[O](27816)'],
    products = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH]=C([CH2])C([CH2])([O])[CH]O(27817)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.854,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cd_H_out_singleH] for rate rule [R5HJ_1;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #6649',
    isomers = [
        '[CH2]C(=C)C([CH2])([O])[CH][O](25001)',
    ],
    reactants = [
        ('C=C([O])[CH][O](2850)', 'allene(458)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #6649',
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


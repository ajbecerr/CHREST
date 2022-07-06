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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07344,0.0423233,-2.94338e-05,1.10648e-08,-1.78046e-12,78231.9,21.8978], Tmin=(100,'K'), Tmax=(1387.7,'K')), NASAPolynomial(coeffs=[8.19337,0.0246827,-1.03656e-05,1.90418e-09,-1.30139e-13,76533.4,-9.63256], Tmin=(1387.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[C]=[CH](745)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.00236431,1.31486e-06,-2.64792e-09,1.00931e-12,101918,6.62296], Tmin=(100,'K'), Tmax=(1016.97,'K')), NASAPolynomial(coeffs=[4.17916,0.00251347,-9.4347e-07,1.68719e-10,-1.1583e-14,101782,4.75425], Tmin=(1016.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C[CH2](21)',
    structure = SMILES('[CH2]C[CH2]'),
    E0 = (289.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00735511,'amu*angstrom^2'), symmetry=1, barrier=(6.35287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00735341,'amu*angstrom^2'), symmetry=1, barrier=(6.35109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2737.71,'J/mol'), sigma=(5.09232,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=427.62 K, Pc=47.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.13019,0.0139359,1.71983e-05,-2.69374e-08,9.94944e-12,34911,13.362], Tmin=(100,'K'), Tmax=(1011.02,'K')), NASAPolynomial(coeffs=[5.48784,0.0173101,-6.65292e-06,1.21645e-09,-8.5034e-14,33785.1,-1.2487], Tmin=(1011.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH]CC[C]=[CH](746)',
    structure = SMILES('[CH]CC[C]=[CH]'),
    E0 = (892.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,3120,650,792.5,1650,180,820.869,1836.35,1841.86],'cm^-1')),
        HinderedRotor(inertia=(0.432526,'amu*angstrom^2'), symmetry=1, barrier=(13.248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.576155,'amu*angstrom^2'), symmetry=1, barrier=(13.2517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0277263,'amu*angstrom^2'), symmetry=1, barrier=(13.2602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98095,0.0446311,-3.83177e-05,1.82196e-08,-3.63961e-12,107458,21.4517], Tmin=(100,'K'), Tmax=(1170.13,'K')), NASAPolynomial(coeffs=[8.57661,0.0220841,-9.41413e-06,1.75195e-09,-1.21224e-13,105914,-11.4047], Tmin=(1170.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(892.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]CC[CH2](747)',
    structure = SMILES('[C]=[C]CC[CH2]'),
    E0 = (960.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,277.933,1701.88],'cm^-1')),
        HinderedRotor(inertia=(0.00218229,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179879,'amu*angstrom^2'), symmetry=1, barrier=(9.86242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179885,'amu*angstrom^2'), symmetry=1, barrier=(9.86243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37388,0.0409338,-2.81886e-05,-5.97432e-09,1.68778e-11,115622,20.8597], Tmin=(100,'K'), Tmax=(530.986,'K')), NASAPolynomial(coeffs=[5.1783,0.0276417,-1.27704e-05,2.4534e-09,-1.72313e-13,115213,8.06531], Tmin=(530.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(960.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
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
    label = '[CH]=CCC=C(332)',
    structure = SMILES('[CH]=CCC=C'),
    E0 = (335.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,249.692],'cm^-1')),
        HinderedRotor(inertia=(0.32851,'amu*angstrom^2'), symmetry=1, barrier=(14.5412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327584,'amu*angstrom^2'), symmetry=1, barrier=(14.5389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15388,0.0316793,6.9434e-06,-2.9126e-08,1.26948e-11,40467.3,19.8537], Tmin=(100,'K'), Tmax=(1001.65,'K')), NASAPolynomial(coeffs=[10.0558,0.0208755,-7.95409e-06,1.47295e-09,-1.04738e-13,37843.3,-23.4778], Tmin=(1001.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CCC(749)',
    structure = SMILES('C#C[CH]CC'),
    E0 = (272.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01894,0.0389033,-1.36053e-05,-1.10514e-08,8.25817e-12,32791.2,17.1163], Tmin=(100,'K'), Tmax=(886.385,'K')), NASAPolynomial(coeffs=[9.03355,0.0209579,-6.437e-06,1.00642e-09,-6.4123e-14,31009.1,-18.9171], Tmin=(886.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH][CH2](487)',
    structure = SMILES('C#CC[CH][CH2]'),
    E0 = (525.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,351.212],'cm^-1')),
        HinderedRotor(inertia=(0.00137256,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00136973,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00131725,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55694,0.0347121,-1.54801e-05,-1.2605e-08,1.47845e-11,63251.6,19.9112], Tmin=(100,'K'), Tmax=(588.934,'K')), NASAPolynomial(coeffs=[5.00988,0.0252601,-9.76533e-06,1.70837e-09,-1.13576e-13,62837.6,8.31453], Tmin=(588.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C=CC[CH2](750)',
    structure = SMILES('C#C[CH]C[CH2]'),
    E0 = (477.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0041749,'amu*angstrom^2'), symmetry=1, barrier=(6.96764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.94423,'amu*angstrom^2'), symmetry=1, barrier=(67.6936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1369,'amu*angstrom^2'), symmetry=1, barrier=(67.7355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76795,0.0430915,-3.64273e-05,1.77222e-08,-3.47094e-12,57486.3,19.7996], Tmin=(100,'K'), Tmax=(1346.56,'K')), NASAPolynomial(coeffs=[9.46879,0.0175599,-5.02767e-06,7.11761e-10,-4.0866e-14,55653.2,-18.7497], Tmin=(1346.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(RCCJ)"""),
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
    label = '[CH]=[C]C[CH][CH2](489)',
    structure = SMILES('[CH]=[C]C[CH][CH2]'),
    E0 = (844.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,1852.76],'cm^-1')),
        HinderedRotor(inertia=(0.101292,'amu*angstrom^2'), symmetry=1, barrier=(2.32891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101805,'amu*angstrom^2'), symmetry=1, barrier=(2.34071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101094,'amu*angstrom^2'), symmetry=1, barrier=(2.32434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20105,0.0433509,-5.18958e-05,4.42318e-08,-1.60035e-11,101610,24.1789], Tmin=(100,'K'), Tmax=(817.043,'K')), NASAPolynomial(coeffs=[3.85249,0.0286472,-1.27502e-05,2.3761e-09,-1.62647e-13,101561,17.8973], Tmin=(817.043,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(844.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C[CH2](752)',
    structure = SMILES('[CH][C]=CC[CH2]'),
    E0 = (761.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,515.587,515.6,515.61,515.611],'cm^-1')),
        HinderedRotor(inertia=(0.280388,'amu*angstrom^2'), symmetry=1, barrier=(52.9053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280431,'amu*angstrom^2'), symmetry=1, barrier=(52.9053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280438,'amu*angstrom^2'), symmetry=1, barrier=(52.9054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0453,0.0405274,-2.17286e-05,5.57859e-09,-5.79171e-13,91626.3,22.277], Tmin=(100,'K'), Tmax=(2104.67,'K')), NASAPolynomial(coeffs=[11.5253,0.0225101,-8.88762e-06,1.51109e-09,-9.60148e-14,87635.9,-30.5132], Tmin=(2104.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(761.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[CH]C(471)',
    structure = SMILES('[CH]=[C]C[CH]C'),
    E0 = (639.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,232.443],'cm^-1')),
        HinderedRotor(inertia=(0.00275455,'amu*angstrom^2'), symmetry=1, barrier=(5.88367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153223,'amu*angstrom^2'), symmetry=1, barrier=(5.88522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153684,'amu*angstrom^2'), symmetry=1, barrier=(5.88253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91189,0.0302195,2.7808e-05,-1.26114e-07,1.15622e-10,76896.3,20.0276], Tmin=(100,'K'), Tmax=(409.644,'K')), NASAPolynomial(coeffs=[3.8483,0.0313221,-1.3748e-05,2.57465e-09,-1.78009e-13,76733.6,15.2963], Tmin=(409.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]C[CH2](490)',
    structure = SMILES('[CH]C=CC[CH2]'),
    E0 = (523.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98638,0.037094,-1.47885e-06,-1.70765e-08,7.2494e-12,63027,22.0103], Tmin=(100,'K'), Tmax=(1082.65,'K')), NASAPolynomial(coeffs=[8.03056,0.0297375,-1.20336e-05,2.19827e-09,-1.51499e-13,60840.7,-11.6824], Tmin=(1082.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]CC(753)',
    structure = SMILES('[CH][C]=CCC'),
    E0 = (555.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85888,0.0414772,-1.68835e-05,2.79592e-10,1.01452e-12,66950.4,21.1656], Tmin=(100,'K'), Tmax=(1338.63,'K')), NASAPolynomial(coeffs=[8.63268,0.0290645,-1.17464e-05,2.08992e-09,-1.39464e-13,64435.5,-16.1094], Tmin=(1338.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC[CH][CH2](330)',
    structure = SMILES('[CH]=CC[CH][CH2]'),
    E0 = (606.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,1500.52],'cm^-1')),
        HinderedRotor(inertia=(0.0957783,'amu*angstrom^2'), symmetry=1, barrier=(2.20213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0966041,'amu*angstrom^2'), symmetry=1, barrier=(2.22112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0972232,'amu*angstrom^2'), symmetry=1, barrier=(2.23535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3052.11,'J/mol'), sigma=(5.53315,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.73 K, Pc=40.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43023,0.0365998,-1.99417e-05,5.4832e-09,-6.43183e-13,72997.5,22.8669], Tmin=(100,'K'), Tmax=(1732.13,'K')), NASAPolynomial(coeffs=[7.07418,0.0258755,-1.06546e-05,1.90871e-09,-1.27271e-13,71388.7,-2.08863], Tmin=(1732.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C[CH][C]=C(279)',
    structure = SMILES('[CH2][C]=CC[CH2]'),
    E0 = (542.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.124939,'amu*angstrom^2'), symmetry=1, barrier=(27.1579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204984,'amu*angstrom^2'), symmetry=1, barrier=(4.71298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00597409,'amu*angstrom^2'), symmetry=1, barrier=(27.137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93243,0.0395004,-1.87977e-05,1.87552e-10,1.75516e-12,65271.7,21.8507], Tmin=(100,'K'), Tmax=(1147.62,'K')), NASAPolynomial(coeffs=[9.44538,0.0226846,-9.06622e-06,1.64916e-09,-1.13131e-13,62930.3,-18.1178], Tmin=(1147.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C[C]=C(278)',
    structure = SMILES('[CH2][CH]C[C]=C'),
    E0 = (597.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,1361.14,1361.33],'cm^-1')),
        HinderedRotor(inertia=(0.157829,'amu*angstrom^2'), symmetry=1, barrier=(3.62879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157827,'amu*angstrom^2'), symmetry=1, barrier=(3.62875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157809,'amu*angstrom^2'), symmetry=1, barrier=(3.62835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86901,0.0335619,-1.56825e-05,3.25552e-09,-2.60722e-13,71861.2,21.3135], Tmin=(100,'K'), Tmax=(2708.28,'K')), NASAPolynomial(coeffs=[12.4387,0.019428,-7.85429e-06,1.32855e-09,-8.28455e-14,66677.7,-34.3892], Tmin=(2708.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S)"""),
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
    E0 = (649.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (895.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1111.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1194.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1104.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1172.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (658.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (713.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (713.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (747.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (700.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (678.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (669.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (928.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1056.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (973.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (801.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (840.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (767.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (791.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (795.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (858.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['C=C(17)', '[CH]=C=C(264)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]C(=C)C[CH2](287)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(15)', '[CH]=[C]C[CH2](308)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=[CH](745)', '[CH2]C[CH2](21)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.45208e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(9)', '[CH]CC[C]=[CH](746)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(9)', '[C]=[C]CC[CH2](747)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=C1CCC1(748)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=CCC=C(332)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=C=CCC(749)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(9)', 'C#CC[CH][CH2](487)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(9)', '[CH]=C=CC[CH2](750)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(17)', '[CH]=[C][CH2](751)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.00310793,'m^3/(mol*s)'), n=2.49201, Ea=(21.4614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH2](18)', '[CH]=C=C(264)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0151074,'m^3/(mol*s)'), n=2.43138, Ea=(27.0165,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH2](18)', '[CH]=[C][CH2](751)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(9)', '[CH]=[C]C[CH][CH2](489)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(9)', '[CH]=[C][CH]C[CH2](752)'],
    products = ['[CH]=[C]CC[CH2](280)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=[C]C[CH]C(471)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=C[CH]C[CH2](490)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.57388e+10,'s^-1'), n=0.905639, Ea=(190.563,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=[C][CH]CC(753)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH]=CC[CH][CH2](330)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH2][CH]C[C]=C(278)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #186',
    isomers = [
        '[CH]=[C]CC[CH2](280)',
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
    label = 'PDepNetwork #186',
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


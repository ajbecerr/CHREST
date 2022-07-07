species(
    label = '[CH2]CC([O])[C]([CH2])[O](25253)',
    structure = SMILES('[CH2]CC([O])[C]([CH2])[O]'),
    E0 = (541.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,485.553,486.159,486.39,3012.19],'cm^-1')),
        HinderedRotor(inertia=(0.0178023,'amu*angstrom^2'), symmetry=1, barrier=(2.98664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193309,'amu*angstrom^2'), symmetry=1, barrier=(32.4199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360371,'amu*angstrom^2'), symmetry=1, barrier=(8.28564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0494057,'amu*angstrom^2'), symmetry=1, barrier=(8.28608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437572,0.0828961,-0.000110456,8.11855e-08,-2.40884e-11,65210.6,30.8209], Tmin=(100,'K'), Tmax=(822.309,'K')), NASAPolynomial(coeffs=[11.3704,0.0297134,-1.34408e-05,2.52997e-09,-1.74646e-13,63412.6,-19.7842], Tmin=(822.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = 'C2H4(22)',
    structure = SMILES('C=C'),
    E0 = (41.9072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98876,-0.00674712,5.04398e-05,-5.70744e-08,2.04945e-11,5047.05,3.80494], Tmin=(100,'K'), Tmax=(946.013,'K')), NASAPolynomial(coeffs=[4.59022,0.00872722,-2.66493e-06,4.8171e-10,-3.6069e-14,4127.02,-3.3246], Tmin=(946.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.9072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]CC([O])C(=C)[O](25257)',
    structure = SMILES('[CH2]CC([O])C(=C)[O]'),
    E0 = (145.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.761561,0.0624558,-4.52804e-05,9.61692e-09,2.02813e-12,17653.3,30.3659], Tmin=(100,'K'), Tmax=(1007.61,'K')), NASAPolynomial(coeffs=[15.8208,0.0197872,-7.23723e-06,1.30222e-09,-9.10164e-14,13749.8,-46.7113], Tmin=(1007.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])([O])[CH][O](24281)',
    structure = SMILES('[CH2]CC([CH2])([O])[CH][O]'),
    E0 = (537.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,292.86,876.722,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151665,'amu*angstrom^2'), symmetry=1, barrier=(3.48708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151665,'amu*angstrom^2'), symmetry=1, barrier=(3.48708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151665,'amu*angstrom^2'), symmetry=1, barrier=(3.48708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151665,'amu*angstrom^2'), symmetry=1, barrier=(3.48708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4235.55,'J/mol'), sigma=(7.22505,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.58 K, Pc=25.48 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0865506,0.097572,-0.000165747,1.50578e-07,-5.24448e-11,64735.7,30.8199], Tmin=(100,'K'), Tmax=(862.188,'K')), NASAPolynomial(coeffs=[7.87782,0.0369948,-1.78531e-05,3.35751e-09,-2.27275e-13,64300.2,-0.347169], Tmin=(862.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(RCCJ)"""),
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
    label = '[CH2]C[CH][C]([CH2])[O](26290)',
    structure = SMILES('[CH2]C[CH][C]([CH2])[O]'),
    E0 = (683.171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,684.289,1828.4,4000],'cm^-1')),
        HinderedRotor(inertia=(0.118275,'amu*angstrom^2'), symmetry=1, barrier=(2.71937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710303,'amu*angstrom^2'), symmetry=1, barrier=(16.3313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00818393,'amu*angstrom^2'), symmetry=1, barrier=(2.71632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77211,'amu*angstrom^2'), symmetry=1, barrier=(40.7443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22813,0.0637334,-7.18828e-05,4.59862e-08,-1.21255e-11,82264,28.3141], Tmin=(100,'K'), Tmax=(913.745,'K')), NASAPolynomial(coeffs=[9.57226,0.0272056,-1.19178e-05,2.23504e-09,-1.54981e-13,80739.1,-11.1887], Tmin=(913.745,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C][O](2821)',
    structure = SMILES('[CH2][C][O]'),
    E0 = (648.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,2091.34],'cm^-1')),
        HinderedRotor(inertia=(0.0328816,'amu*angstrom^2'), symmetry=1, barrier=(10.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26468,0.0187769,-3.12134e-05,3.2382e-08,-1.33361e-11,78037.6,11.2589], Tmin=(100,'K'), Tmax=(748.797,'K')), NASAPolynomial(coeffs=[3.86755,0.0107985,-5.69995e-06,1.18131e-09,-8.60475e-14,78080.7,9.41548], Tmin=(748.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C[CH][O](1549)',
    structure = SMILES('[CH2]C[CH][O]'),
    E0 = (338.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1524.15,1525.95],'cm^-1')),
        HinderedRotor(inertia=(0.310045,'amu*angstrom^2'), symmetry=1, barrier=(7.12854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00431249,'amu*angstrom^2'), symmetry=1, barrier=(7.11449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81289,0.029873,-3.56716e-05,3.591e-08,-1.53006e-11,40735.5,17.4078], Tmin=(100,'K'), Tmax=(773.473,'K')), NASAPolynomial(coeffs=[2.08215,0.0257045,-1.21752e-05,2.37389e-09,-1.67404e-13,41086.3,22.2823], Tmin=(773.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(RCCJ) + radical(CCsJOH)"""),
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
    label = '[CH2][C]([O])C([CH2])[O](1982)',
    structure = SMILES('[CH2][C]([O])C([CH2])[O]'),
    E0 = (571.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,380.245,722.45,2796.96],'cm^-1')),
        HinderedRotor(inertia=(0.481164,'amu*angstrom^2'), symmetry=1, barrier=(11.0629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.479644,'amu*angstrom^2'), symmetry=1, barrier=(11.028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768815,'amu*angstrom^2'), symmetry=1, barrier=(28.1179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827831,0.073441,-0.000109338,8.49615e-08,-2.59656e-11,68820.2,26.0158], Tmin=(100,'K'), Tmax=(832.925,'K')), NASAPolynomial(coeffs=[11.668,0.0194549,-8.64352e-06,1.58734e-09,-1.06978e-13,67081.3,-23.8985], Tmin=(832.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(CJCO)"""),
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
    label = '[CH]CC([O])[C]([CH2])[O](26291)',
    structure = SMILES('[CH]CC([O])[C]([CH2])[O]'),
    E0 = (784.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301748,0.0857996,-0.000121729,9.17199e-08,-2.74631e-11,94437.9,30.5234], Tmin=(100,'K'), Tmax=(819.5,'K')), NASAPolynomial(coeffs=[12.6218,0.0256666,-1.16647e-05,2.18461e-09,-1.49809e-13,92418.6,-26.4615], Tmin=(819.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(784.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]([O])C([O])C[CH2](26292)',
    structure = SMILES('[CH][C]([O])C([O])C[CH2]'),
    E0 = (777.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528843,0.0808723,-0.000110647,8.23151e-08,-2.45982e-11,93666.9,30.863], Tmin=(100,'K'), Tmax=(817.856,'K')), NASAPolynomial(coeffs=[11.4636,0.0273956,-1.25739e-05,2.37637e-09,-1.64233e-13,91878.2,-19.6928], Tmin=(817.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(777.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]CC1OC1([CH2])[O](26293)',
    structure = SMILES('[CH2]CC1OC1([CH2])[O]'),
    E0 = (277.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0787014,0.0724913,-7.04961e-05,3.60293e-08,-7.03879e-12,33566.7,27.4328], Tmin=(100,'K'), Tmax=(1433.94,'K')), NASAPolynomial(coeffs=[16.1793,0.0180897,-3.66233e-06,3.42129e-10,-1.2351e-14,29924.8,-52.6445], Tmin=(1433.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CJC(O)2C) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([O])C1CCO1(25261)',
    structure = SMILES('[CH2][C]([O])C1CCO1'),
    E0 = (289.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441924,0.0664436,-6.17088e-05,3.13478e-08,-6.16225e-12,34967.2,26.606], Tmin=(100,'K'), Tmax=(1422.18,'K')), NASAPolynomial(coeffs=[13.8222,0.0205566,-4.60548e-06,4.98996e-10,-2.20935e-14,31996.1,-39.7238], Tmin=(1422.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1([O])CCC1[O](26281)',
    structure = SMILES('[CH2]C1([O])CCC1[O]'),
    E0 = (298.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897503,0.0560152,-1.95926e-05,-1.82798e-08,1.18126e-11,36019.1,24.0508], Tmin=(100,'K'), Tmax=(993.842,'K')), NASAPolynomial(coeffs=[16.4027,0.0209105,-7.81382e-06,1.459e-09,-1.05454e-13,31588.9,-57.4404], Tmin=(993.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]CC(O)=C([CH2])[O](26294)',
    structure = SMILES('[CH2]CC(O)=C([CH2])[O]'),
    E0 = (20.3033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.289795,0.0783747,-8.31499e-05,4.29892e-08,-8.33427e-12,2610.06,31.3737], Tmin=(100,'K'), Tmax=(1438.47,'K')), NASAPolynomial(coeffs=[20.5916,0.00999387,-1.08761e-06,-2.81451e-11,8.35349e-15,-2330.15,-73.2493], Tmin=(1438.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.3033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])[O](26295)',
    structure = SMILES('[CH2]CC(=O)C([CH2])[O]'),
    E0 = (213.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944442,0.0706604,-8.00766e-05,5.09559e-08,-1.33606e-11,25775.4,30.1543], Tmin=(100,'K'), Tmax=(917.941,'K')), NASAPolynomial(coeffs=[10.3042,0.029875,-1.34303e-05,2.55389e-09,-1.78596e-13,24057.1,-14.2], Tmin=(917.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C=OCOJ) + radical(CJCC=O) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])C(=O)CC(26296)',
    structure = SMILES('[CH2]C([O])=C([O])CC'),
    E0 = (-47.1383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247926,0.0712258,-7.17517e-05,3.69308e-08,-7.30628e-12,-5524.7,29.157], Tmin=(100,'K'), Tmax=(1354,'K')), NASAPolynomial(coeffs=[17.319,0.0153303,-3.77616e-06,4.81527e-10,-2.60759e-14,-9646.7,-56.5252], Tmin=(1354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.1383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][C]([O])C(O)C=C(26297)',
    structure = SMILES('[CH2][C]([O])C(O)C=C'),
    E0 = (232.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617169,0.0749332,-8.33472e-05,4.87447e-08,-1.14088e-11,28144.7,29.9615], Tmin=(100,'K'), Tmax=(1037.1,'K')), NASAPolynomial(coeffs=[13.6556,0.0246458,-1.06152e-05,1.99179e-09,-1.38798e-13,25440.2,-33.4164], Tmin=(1037.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])C([O])C=C(26298)',
    structure = SMILES('[CH2]C([O])C([O])C=C'),
    E0 = (286.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582469,0.0663766,-5.50375e-05,2.04652e-08,-2.25608e-12,34616.4,30.4667], Tmin=(100,'K'), Tmax=(1099.34,'K')), NASAPolynomial(coeffs=[16.3838,0.0203947,-8.0051e-06,1.46921e-09,-1.02434e-13,30446.5,-50.426], Tmin=(1099.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]CC(=O)[C]([CH2])[O](26299)',
    structure = SMILES('[CH2]CC([O])=C([CH2])[O]'),
    E0 = (158.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,440.78,440.792,3966.86],'cm^-1')),
        HinderedRotor(inertia=(0.000867405,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101974,'amu*angstrom^2'), symmetry=1, barrier=(14.0699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102058,'amu*angstrom^2'), symmetry=1, barrier=(14.0692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424313,0.0703666,-7.68529e-05,4.25472e-08,-9.0392e-12,19151.9,30.3074], Tmin=(100,'K'), Tmax=(1243.89,'K')), NASAPolynomial(coeffs=[16.7371,0.0136697,-3.37005e-06,4.23839e-10,-2.2464e-14,15421.6,-50.6343], Tmin=(1243.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][C]([O])C([O])C=C(25467)',
    structure = SMILES('[CH2][C]([O])C([O])C=C'),
    E0 = (463.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,434.42,434.421,434.43,434.452],'cm^-1')),
        HinderedRotor(inertia=(0.000893192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092955,'amu*angstrom^2'), symmetry=1, barrier=(12.4524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0929712,'amu*angstrom^2'), symmetry=1, barrier=(12.4522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773567,0.0710285,-7.78057e-05,4.42375e-08,-1.00443e-11,55845.5,29.2072], Tmin=(100,'K'), Tmax=(1068.26,'K')), NASAPolynomial(coeffs=[13.6818,0.0226944,-9.9365e-06,1.8821e-09,-1.31955e-13,53087.6,-33.9196], Tmin=(1068.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = 'CH2CH2CHO(1580)',
    structure = SMILES('[CH2]CC=O'),
    E0 = (11.2619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.221236,'amu*angstrom^2'), symmetry=1, barrier=(5.08666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221178,'amu*angstrom^2'), symmetry=1, barrier=(5.08531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76345,0.0293577,-2.47893e-05,1.49241e-08,-4.38506e-12,1397.08,14.4322], Tmin=(100,'K'), Tmax=(767.822,'K')), NASAPolynomial(coeffs=[4.24224,0.0216538,-9.7387e-06,1.85597e-09,-1.3004e-13,1170,7.68863], Tmin=(767.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.2619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH2CH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[CH]C(=C)[O](2640)',
    structure = SMILES('[CH2]CC=C([CH2])[O]'),
    E0 = (234.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,393.035,400.485],'cm^-1')),
        HinderedRotor(inertia=(0.0011284,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.044464,'amu*angstrom^2'), symmetry=1, barrier=(5.0623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0462043,'amu*angstrom^2'), symmetry=1, barrier=(5.19959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32742,0.0524088,-3.53875e-05,5.62558e-09,2.75441e-12,28349.3,25.5104], Tmin=(100,'K'), Tmax=(962.505,'K')), NASAPolynomial(coeffs=[12.62,0.0197354,-6.6865e-06,1.13536e-09,-7.64686e-14,25515.1,-31.9689], Tmin=(962.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C2H4(T)(23)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,995.342,995.492,2246.13,2948.67],'cm^-1')),
        HinderedRotor(inertia=(0.00723049,'amu*angstrom^2'), symmetry=1, barrier=(5.08405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40735,0.0100313,6.40897e-06,-1.41287e-08,5.92653e-12,38288.2,6.11706], Tmin=(100,'K'), Tmax=(954.266,'K')), NASAPolynomial(coeffs=[5.52251,0.00856169,-2.90741e-06,5.02347e-10,-3.44567e-14,37547.7,-5.75289], Tmin=(954.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C[C]([O])[C]([CH2])[O](26300)',
    structure = SMILES('[CH2]C[C]([O])[C]([CH2])[O]'),
    E0 = (717.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,363.333,366.667,370,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,289.279,290.435,1229.33,1229.95],'cm^-1')),
        HinderedRotor(inertia=(0.00199293,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118857,'amu*angstrom^2'), symmetry=1, barrier=(7.16639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12053,'amu*angstrom^2'), symmetry=1, barrier=(7.16353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12026,'amu*angstrom^2'), symmetry=1, barrier=(7.17046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145497,0.0933952,-0.000154361,1.33335e-07,-4.44157e-11,86460.5,31.2845], Tmin=(100,'K'), Tmax=(866.099,'K')), NASAPolynomial(coeffs=[10.436,0.0290054,-1.36368e-05,2.53271e-09,-1.70134e-13,85310.5,-13.2302], Tmin=(866.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]C([O])[C]([CH2])[O](26301)',
    structure = SMILES('[CH2][CH]C([O])[C]([CH2])[O]'),
    E0 = (741.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,2627.13,2627.49],'cm^-1')),
        HinderedRotor(inertia=(0.158035,'amu*angstrom^2'), symmetry=1, barrier=(3.63354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0291949,'amu*angstrom^2'), symmetry=1, barrier=(16.048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292484,'amu*angstrom^2'), symmetry=1, barrier=(16.0496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0606058,'amu*angstrom^2'), symmetry=1, barrier=(33.2552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627885,0.0777299,-0.000101971,7.17972e-08,-2.02272e-11,89247.4,32.6344], Tmin=(100,'K'), Tmax=(867.41,'K')), NASAPolynomial(coeffs=[11.9953,0.0253096,-1.13208e-05,2.12588e-09,-1.46858e-13,87275.4,-20.5899], Tmin=(867.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCJCO) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C](O)[C]([CH2])[O](26302)',
    structure = SMILES('[CH2]C[C](O)[C]([CH2])[O]'),
    E0 = (487.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0168361,0.0973798,-0.000160211,1.38247e-07,-4.59385e-11,58759.9,32.0592], Tmin=(100,'K'), Tmax=(872.214,'K')), NASAPolynomial(coeffs=[10.5237,0.030759,-1.41996e-05,2.6147e-09,-1.7466e-13,57616.6,-13.3655], Tmin=(872.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C[C]([O])C([CH2])[O](26303)',
    structure = SMILES('[CH2]C[C]([O])C([CH2])[O]'),
    E0 = (541.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,402.931,2304.26],'cm^-1')),
        HinderedRotor(inertia=(0.00332657,'amu*angstrom^2'), symmetry=1, barrier=(12.5486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545248,'amu*angstrom^2'), symmetry=1, barrier=(12.5363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0240754,'amu*angstrom^2'), symmetry=1, barrier=(12.5483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0570029,'amu*angstrom^2'), symmetry=1, barrier=(29.6856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437572,0.0828961,-0.000110456,8.11855e-08,-2.40884e-11,65210.6,30.8209], Tmin=(100,'K'), Tmax=(822.309,'K')), NASAPolynomial(coeffs=[11.3704,0.0297134,-1.34408e-05,2.52997e-09,-1.74646e-13,63412.6,-19.7842], Tmin=(822.309,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])C([O])[CH]C(26304)',
    structure = SMILES('[CH2][C]([O])C([O])[CH]C'),
    E0 = (535.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640362,0.0764061,-8.9428e-05,5.67653e-08,-1.45602e-11,64562.6,30.8034], Tmin=(100,'K'), Tmax=(945.846,'K')), NASAPolynomial(coeffs=[12.1686,0.0276544,-1.21157e-05,2.27426e-09,-1.57913e-13,62381.7,-24.1721], Tmin=(945.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CCJCO) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]C(O)[C]([CH2])[O](26305)',
    structure = SMILES('[CH2][CH]C(O)[C]([CH2])[O]'),
    E0 = (510.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439403,0.0820448,-0.000109088,7.85344e-08,-2.26206e-11,61548,33.5015], Tmin=(100,'K'), Tmax=(849.33,'K')), NASAPolynomial(coeffs=[12.1322,0.0269767,-1.18325e-05,2.19559e-09,-1.50351e-13,59561.8,-20.9999], Tmin=(849.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]C([O])C([CH2])[O](26306)',
    structure = SMILES('[CH2][CH]C([O])C([CH2])[O]'),
    E0 = (564.434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,814.395,3643.5],'cm^-1')),
        HinderedRotor(inertia=(1.09373,'amu*angstrom^2'), symmetry=1, barrier=(25.147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0534382,'amu*angstrom^2'), symmetry=1, barrier=(25.146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00617712,'amu*angstrom^2'), symmetry=1, barrier=(2.90683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0061777,'amu*angstrom^2'), symmetry=1, barrier=(2.90669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6883,0.0702034,-6.9542e-05,3.59563e-08,-7.44455e-12,68007.2,32.9853], Tmin=(100,'K'), Tmax=(1166.49,'K')), NASAPolynomial(coeffs=[14.2415,0.0237279,-9.77861e-06,1.8005e-09,-1.24312e-13,64845.3,-34.4885], Tmin=(1166.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CCJCO) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([O])[C]([O])CC(26307)',
    structure = SMILES('[CH2][C]([O])[C]([O])CC'),
    E0 = (512.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.12474,0.0924608,-0.000143185,1.20136e-07,-3.95853e-11,61777.1,29.573], Tmin=(100,'K'), Tmax=(849.215,'K')), NASAPolynomial(coeffs=[10.5716,0.0314241,-1.44787e-05,2.69302e-09,-1.82232e-13,60429.3,-16.6081], Tmin=(849.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C[C]([O])[C]([CH2])O(26308)',
    structure = SMILES('[CH2]C[C]([O])[C]([CH2])O'),
    E0 = (487.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0168361,0.0973798,-0.000160211,1.38247e-07,-4.59385e-11,58759.9,32.0592], Tmin=(100,'K'), Tmax=(872.214,'K')), NASAPolynomial(coeffs=[10.5237,0.030759,-1.41996e-05,2.6147e-09,-1.7466e-13,57616.6,-13.3655], Tmin=(872.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C[C]([O])[C](C)[O](26309)',
    structure = SMILES('[CH2]C[C]([O])[C](C)[O]'),
    E0 = (506.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.328261,0.087842,-0.00013334,1.12589e-07,-3.76427e-11,61007,29.9953], Tmin=(100,'K'), Tmax=(838.833,'K')), NASAPolynomial(coeffs=[9.46672,0.0330565,-1.53297e-05,2.8706e-09,-1.95452e-13,59868.2,-10.136], Tmin=(838.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C2CsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([O])[C]([CH2])O(26310)',
    structure = SMILES('[CH2][CH]C([O])[C]([CH2])O'),
    E0 = (510.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439405,0.0820448,-0.000109088,7.85343e-08,-2.26206e-11,61548,33.5015], Tmin=(100,'K'), Tmax=(849.342,'K')), NASAPolynomial(coeffs=[12.1322,0.0269767,-1.18325e-05,2.19559e-09,-1.50351e-13,59561.8,-20.9999], Tmin=(849.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCJCO) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([O])[C](C)[O](26311)',
    structure = SMILES('[CH2][CH]C([O])[C](C)[O]'),
    E0 = (529.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.869053,0.0714804,-7.84463e-05,4.76106e-08,-1.18497e-11,63791.4,31.1359], Tmin=(100,'K'), Tmax=(966.744,'K')), NASAPolynomial(coeffs=[11.0698,0.0292735,-1.29577e-05,2.44948e-09,-1.7092e-13,61819.1,-17.7317], Tmin=(966.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCJCO) + radical(RCCJ)"""),
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
    E0 = (541.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (541.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (704.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1202.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1042.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1009.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (995.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (989.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (544.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (549.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (549.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (564.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (564.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (604.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (604.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (604.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (541.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (685.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (715.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (541.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (541.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (635.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (889.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (929.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (952.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (655.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (699.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (693.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (616.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (680.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (659.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (651.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (698.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (664.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (672.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['C2H4(22)', 'C=C([O])[CH][O](2850)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]CC([O])C(=C)[O](25257)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]CC([CH2])([O])[CH][O](24281)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH2]C[CH][C]([CH2])[O](26290)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][C][O](2821)', '[CH2]C[CH][O](1549)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][C]([O])C([CH2])[O](1982)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH]CC([O])[C]([CH2])[O](26291)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH][C]([O])C([O])C[CH2](26292)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]CC1OC1([CH2])[O](26293)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][C]([O])C1CCO1(25261)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]C1([O])CCC1[O](26281)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]CC(O)=C([CH2])[O](26294)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]CC(=O)C([CH2])[O](26295)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][C]([O])C(=O)CC(26296)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][C]([O])C(O)C=C(26297)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]C([O])C([O])C=C(26298)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]CC(=O)[C]([CH2])[O](26299)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(171.259,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 167.7 to 171.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][C]([O])C([O])C=C(25467)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][C][O](2821)', 'CH2CH2CHO(1580)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O(4)', '[CH2]C[CH]C(=C)[O](2640)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(63.2931,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 63.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C2H4(T)(23)', 'C=C([O])[CH][O](2850)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4679.9,'m^3/(mol*s)'), n=0.573452, Ea=(202.057,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 198.3 to 202.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C2H4(22)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00304282,'m^3/(mol*s)'), n=2.49542, Ea=(21.4551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C2H4(T)(23)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.23391e+07,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C[C]([O])[C]([CH2])[O](26300)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH]C([O])[C]([CH2])[O](26301)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]C[C](O)[C]([CH2])[O](26302)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C[C]([O])C([CH2])[O](26303)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][C]([O])C([O])[CH]C(26304)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][CH]C(O)[C]([CH2])[O](26305)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]C([O])C([CH2])[O](26306)'],
    products = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][C]([O])[C]([O])CC(26307)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]C[C]([O])[C]([CH2])O(26308)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.23689e+09,'s^-1'), n=1.09705, Ea=(110.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;XH_out] for rate rule [R3HJ;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2]C[C]([O])[C](C)[O](26309)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][CH]C([O])[C]([CH2])O(26310)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.04154e+06,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]CC([O])[C]([CH2])[O](25253)'],
    products = ['[CH2][CH]C([O])[C](C)[O](26311)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #5843',
    isomers = [
        '[CH2]CC([O])[C]([CH2])[O](25253)',
    ],
    reactants = [
        ('C2H4(22)', 'C=C([O])[CH][O](2850)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #5843',
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


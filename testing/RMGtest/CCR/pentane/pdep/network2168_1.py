species(
    label = '[CH2][CH]C([CH2])CC=C[CH2](5800)',
    structure = SMILES('[CH2][CH]C([CH2])CC=C[CH2]'),
    E0 = (626.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117574,0.0787344,-5.46928e-05,2.01205e-08,-3.07431e-12,75530.1,38.683], Tmin=(100,'K'), Tmax=(1501.78,'K')), NASAPolynomial(coeffs=[14.8195,0.0395759,-1.55808e-05,2.75806e-09,-1.84022e-13,71114.3,-38.2241], Tmin=(1501.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C=C(3743)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210055,'amu*angstrom^2'), symmetry=1, barrier=(25.2323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779472,'amu*angstrom^2'), symmetry=1, barrier=(93.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56316,0.022343,1.87062e-05,-3.93092e-08,1.63979e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.267,'K')), NASAPolynomial(coeffs=[9.83,0.0151965,-5.22268e-06,9.67646e-10,-7.07852e-14,30607.7,-26.9852], Tmin=(974.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'butadiene13(1350)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30711,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80606,0.0102576,6.17291e-05,-9.01684e-08,3.59136e-11,11658.5,12.0619], Tmin=(100,'K'), Tmax=(946.037,'K')), NASAPolynomial(coeffs=[12.4692,0.0100558,-2.41229e-06,4.5713e-10,-3.93205e-14,8010.87,-43.6362], Tmin=(946.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C=C)C[CH]C=C(4934)',
    structure = SMILES('[CH2]C=CCC([CH2])C=C'),
    E0 = (352.504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3602.68,'J/mol'), sigma=(6.38391,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.73 K, Pc=31.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16257,0.0716653,-2.41815e-05,-1.78449e-08,1.16684e-11,42545.7,34.6823], Tmin=(100,'K'), Tmax=(1003.19,'K')), NASAPolynomial(coeffs=[15.6153,0.0374088,-1.38662e-05,2.48385e-09,-1.72003e-13,38068.7,-46.7784], Tmin=(1003.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]CCC=C[CH2](5895)',
    structure = SMILES('[CH2][CH][CH]CCC=C[CH2]'),
    E0 = (621.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,340.996,979.256,1296.19,2300.1],'cm^-1')),
        HinderedRotor(inertia=(0.0821972,'amu*angstrom^2'), symmetry=1, barrier=(2.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0821972,'amu*angstrom^2'), symmetry=1, barrier=(2.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0821972,'amu*angstrom^2'), symmetry=1, barrier=(2.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0821972,'amu*angstrom^2'), symmetry=1, barrier=(2.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0821972,'amu*angstrom^2'), symmetry=1, barrier=(2.8753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0821972,'amu*angstrom^2'), symmetry=1, barrier=(2.8753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3576.97,'J/mol'), sigma=(6.55299,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.71 K, Pc=28.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66687,0.0617671,1.2969e-05,-1.28943e-07,1.16226e-10,74853.9,34.5884], Tmin=(100,'K'), Tmax=(450.163,'K')), NASAPolynomial(coeffs=[3.81145,0.0581095,-2.61528e-05,4.98082e-09,-3.48892e-13,74504.8,24.2208], Tmin=(450.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH2](502)',
    structure = SMILES('[CH][CH2]'),
    E0 = (557.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1433.18,1433.5],'cm^-1')),
        HinderedRotor(inertia=(0.00559429,'amu*angstrom^2'), symmetry=1, barrier=(8.15686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7298,0.00311091,1.59755e-05,-2.26123e-08,9.32732e-12,67119.1,8.34543], Tmin=(100,'K'), Tmax=(875.199,'K')), NASAPolynomial(coeffs=[4.97018,0.00511238,-6.01253e-07,2.8794e-11,-6.01916e-16,66608.3,0.848354], Tmin=(875.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]CC=C[CH2](652)',
    structure = SMILES('[CH2][CH]CC=C[CH2]'),
    E0 = (474.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,348.938,2083.29],'cm^-1')),
        HinderedRotor(inertia=(0.107987,'amu*angstrom^2'), symmetry=1, barrier=(9.35741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0235326,'amu*angstrom^2'), symmetry=1, barrier=(56.8317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108208,'amu*angstrom^2'), symmetry=1, barrier=(9.35905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0127607,'amu*angstrom^2'), symmetry=1, barrier=(30.8466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44643,0.0505193,-2.85643e-05,7.98419e-09,-9.07072e-13,57209.5,28.1169], Tmin=(100,'K'), Tmax=(1965.55,'K')), NASAPolynomial(coeffs=[13.3496,0.0262965,-1.00793e-05,1.71471e-09,-1.09676e-13,52530.1,-37.3531], Tmin=(1965.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]CC=C[CH2](5280)',
    structure = SMILES('[CH2][CH][CH]CC=C[CH2]'),
    E0 = (645.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,239.279,1036.11,2991.11],'cm^-1')),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36806,0.0623164,-5.56022e-05,4.06775e-08,-1.47416e-11,77731.4,33.1934], Tmin=(100,'K'), Tmax=(737.489,'K')), NASAPolynomial(coeffs=[3.19352,0.0490934,-2.19507e-05,4.14967e-09,-2.88591e-13,77552.5,25.555], Tmin=(737.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH][CH]CC([CH2])C=C(5247)',
    structure = SMILES('[CH][CH]CC([CH2])C=C'),
    E0 = (751.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884291,0.0668622,-5.10869e-05,2.24695e-08,-4.23142e-12,90547.3,33.2699], Tmin=(100,'K'), Tmax=(1225.24,'K')), NASAPolynomial(coeffs=[10.0292,0.0370073,-1.45369e-05,2.58226e-09,-1.73588e-13,88306.3,-12.7065], Tmin=(1225.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][C]C([CH2])CC=C[CH2](9211)',
    structure = SMILES('[CH2][C]C([CH2])CC=C[CH2]'),
    E0 = (880.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27102,0.084735,-6.8987e-05,2.9918e-08,-5.26353e-12,106056,37.427], Tmin=(100,'K'), Tmax=(1354,'K')), NASAPolynomial(coeffs=[16.9296,0.0339209,-1.26941e-05,2.20131e-09,-1.45999e-13,101398,-50.7691], Tmin=(1354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(880.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])CC=C[CH2](9212)',
    structure = SMILES('[CH]C([CH][CH2])CC=C[CH2]'),
    E0 = (869.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16045,0.0787782,-5.70852e-05,2.12929e-08,-3.26374e-12,104770,37.7507], Tmin=(100,'K'), Tmax=(1504.06,'K')), NASAPolynomial(coeffs=[15.8897,0.0369468,-1.53668e-05,2.80148e-09,-1.90171e-13,100038,-44.5542], Tmin=(1504.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(869.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=CCC([CH2])[CH][CH2](9213)',
    structure = SMILES('[CH]C=CCC([CH2])[CH][CH2]'),
    E0 = (845.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450288,0.0774094,-5.05099e-05,1.77937e-08,-2.70986e-12,101875,38.3047], Tmin=(100,'K'), Tmax=(1439.83,'K')), NASAPolynomial(coeffs=[11.0589,0.0479373,-1.9806e-05,3.57715e-09,-2.41401e-13,98819.9,-16.7428], Tmin=(1439.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])CC=C[CH2](9214)',
    structure = SMILES('[CH][CH]C([CH2])CC=C[CH2]'),
    E0 = (869.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0894073,0.0803928,-6.17535e-05,2.54368e-08,-4.33215e-12,104753,37.9981], Tmin=(100,'K'), Tmax=(1370.6,'K')), NASAPolynomial(coeffs=[14.7408,0.0376335,-1.49573e-05,2.6749e-09,-1.80326e-13,100737,-37.3054], Tmin=(1370.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(869.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC1CC1[CH2](5686)',
    structure = SMILES('[CH2]C=CCC1CC1[CH2]'),
    E0 = (377.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576256,0.0550522,3.25012e-05,-8.27968e-08,3.62504e-11,45509.1,31.2689], Tmin=(100,'K'), Tmax=(956.765,'K')), NASAPolynomial(coeffs=[17.8511,0.0327555,-1.08146e-05,1.92478e-09,-1.38164e-13,39918.4,-63.2506], Tmin=(956.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1CC=CCC1[CH2](5688)',
    structure = SMILES('[CH2]C1CC=CCC1[CH2]'),
    E0 = (323.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07509,0.042526,6.40658e-05,-1.12015e-07,4.60377e-11,39003.4,27.1208], Tmin=(100,'K'), Tmax=(947.233,'K')), NASAPolynomial(coeffs=[15.7381,0.0353326,-1.12051e-05,1.95421e-09,-1.39498e-13,33770.4,-55.7845], Tmin=(947.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1C[CH][CH]CC1(5712)',
    structure = SMILES('C=CC1C[CH][CH]CC1'),
    E0 = (294.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59866,0.0397858,4.31611e-05,-6.28985e-08,2.07695e-11,35530.9,27.4109], Tmin=(100,'K'), Tmax=(1128.45,'K')), NASAPolynomial(coeffs=[7.16634,0.0553007,-2.4319e-05,4.64942e-09,-3.27992e-13,32029.9,-10.0675], Tmin=(1128.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(cyclohexane) + radical(cyclohexane)"""),
)

species(
    label = '[CH2]C=CCC(=C)C[CH2](9215)',
    structure = SMILES('[CH2]C=CCC(=C)C[CH2]'),
    E0 = (347.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0383029,0.0748076,-3.48372e-05,-4.44043e-09,5.9914e-12,41980.1,34.5898], Tmin=(100,'K'), Tmax=(1076.25,'K')), NASAPolynomial(coeffs=[15.9134,0.0382686,-1.52185e-05,2.79938e-09,-1.94907e-13,37262,-49.2089], Tmin=(1076.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C](C)C=C(5697)',
    structure = SMILES('[CH2]C=CC[C](C)C=C'),
    E0 = (279.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418244,0.0646102,-6.20513e-06,-3.32799e-08,1.60783e-11,33745.6,30.2376], Tmin=(100,'K'), Tmax=(1017.34,'K')), NASAPolynomial(coeffs=[14.7322,0.0395365,-1.52465e-05,2.79589e-09,-1.96104e-13,29218.3,-47.0017], Tmin=(1017.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)CC=CC(5706)',
    structure = SMILES('[CH2]C=C([CH2])CC=CC'),
    E0 = (280.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,815.112],'cm^-1')),
        HinderedRotor(inertia=(0.92666,'amu*angstrom^2'), symmetry=1, barrier=(21.3057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2885,'amu*angstrom^2'), symmetry=1, barrier=(6.63318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928321,'amu*angstrom^2'), symmetry=1, barrier=(21.3439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926574,'amu*angstrom^2'), symmetry=1, barrier=(21.3038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0102597,'amu*angstrom^2'), symmetry=1, barrier=(21.3315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00300529,0.0752869,-3.39848e-05,-6.02929e-09,6.65549e-12,33908.1,32.1039], Tmin=(100,'K'), Tmax=(1070.02,'K')), NASAPolynomial(coeffs=[16.1077,0.0385292,-1.53232e-05,2.82116e-09,-1.96673e-13,29119.4,-52.9535], Tmin=(1070.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])CC=C=C(9014)',
    structure = SMILES('[CH2]CC([CH2])CC=C=C'),
    E0 = (457.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0367449,0.0821723,-6.14308e-05,2.48071e-08,-4.15126e-12,55157.2,35.4513], Tmin=(100,'K'), Tmax=(1392.3,'K')), NASAPolynomial(coeffs=[14.8875,0.039296,-1.5238e-05,2.68891e-09,-1.79757e-13,51001.4,-41.4889], Tmin=(1392.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]=CCC(C)C=C(5705)',
    structure = SMILES('[CH2][C]=CCC(C)C=C'),
    E0 = (385.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0488069,0.0752304,-3.82878e-05,5.72882e-10,3.94223e-12,46488.5,33.7882], Tmin=(100,'K'), Tmax=(1104.32,'K')), NASAPolynomial(coeffs=[15.4611,0.038935,-1.55155e-05,2.84008e-09,-1.96477e-13,41893.6,-47.4885], Tmin=(1104.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])[CH][CH2](9216)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])[CH][CH2]'),
    E0 = (957.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,727.223,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162386,'amu*angstrom^2'), symmetry=1, barrier=(13.9515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437034,0.0857318,-0.000100769,8.55397e-08,-3.11901e-11,115262,44.2945], Tmin=(100,'K'), Tmax=(805.211,'K')), NASAPolynomial(coeffs=[3.64616,0.0576663,-2.59021e-05,4.85545e-09,-3.33828e-13,115138,31.9481], Tmin=(805.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(957.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C]([CH2])[CH][CH2](9217)',
    structure = SMILES('[CH2][CH]CC[C]([CH2])[CH][CH2]'),
    E0 = (948.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2039.07,2563,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0251275,'amu*angstrom^2'), symmetry=1, barrier=(8.46256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275705,0.0950711,-0.000138818,1.32736e-07,-4.9721e-11,114165,43.4495], Tmin=(100,'K'), Tmax=(855.402,'K')), NASAPolynomial(coeffs=[0.753533,0.0624995,-2.85036e-05,5.30043e-09,-3.59331e-13,115193,47.7066], Tmin=(855.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(948.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C[CH]C[CH2](9218)',
    structure = SMILES('[CH2][CH][C]([CH2])C[CH]C[CH2]'),
    E0 = (948.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2037.98,2613.96,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223837,'amu*angstrom^2'), symmetry=1, barrier=(8.38701,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.316759,0.0955002,-0.000143839,1.40887e-07,-5.34588e-11,114164,43.1304], Tmin=(100,'K'), Tmax=(856.437,'K')), NASAPolynomial(coeffs=[-0.338416,0.0644498,-2.97129e-05,5.5445e-09,-3.76261e-13,115527,53.4931], Tmin=(856.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(948.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])CC1[CH]C1(9219)',
    structure = SMILES('[CH2][CH]C([CH2])CC1[CH]C1'),
    E0 = (742.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507441,0.066368,-2.46e-05,-6.90934e-09,5.29154e-12,89444.6,38.5894], Tmin=(100,'K'), Tmax=(1117.72,'K')), NASAPolynomial(coeffs=[12.3768,0.0421894,-1.67085e-05,3.0308e-09,-2.07892e-13,85648.3,-25.1072], Tmin=(1117.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(742.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(cyclopropane) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])C1[CH2](9206)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])C1[CH2]'),
    E0 = (619.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15066,0.0438318,5.40153e-05,-1.01141e-07,4.29859e-11,74646.6,34.2454], Tmin=(100,'K'), Tmax=(924.085,'K')), NASAPolynomial(coeffs=[14.3444,0.0350063,-1.0036e-05,1.61198e-09,-1.09966e-13,70146.6,-39.5197], Tmin=(924.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(619.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1C[CH]C([CH2])C1(9220)',
    structure = SMILES('[CH2][CH]C1C[CH]C([CH2])C1'),
    E0 = (618.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26418,0.0432938,4.41736e-05,-7.92225e-08,3.107e-11,74479.6,35.5898], Tmin=(100,'K'), Tmax=(990.487,'K')), NASAPolynomial(coeffs=[11.7822,0.0418268,-1.57092e-05,2.88343e-09,-2.04133e-13,70384.4,-25.2078], Tmin=(990.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1[CH]CC([CH2])[CH]C1(5663)',
    structure = SMILES('[CH2]C1[CH]CC([CH2])[CH]C1'),
    E0 = (591.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19343,0.0457444,3.97712e-05,-6.96517e-08,2.56963e-11,71284.3,30.2264], Tmin=(100,'K'), Tmax=(1046.66,'K')), NASAPolynomial(coeffs=[10.4342,0.0490322,-2.02643e-05,3.82612e-09,-2.71211e-13,67235.5,-24.8776], Tmin=(1046.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(591.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclohexane) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])C1[CH2](9209)',
    structure = SMILES('[CH2][CH]C1CC([CH2])C1[CH2]'),
    E0 = (701.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789008,0.0545159,2.19334e-05,-6.38789e-08,2.78379e-11,84510.7,35.3546], Tmin=(100,'K'), Tmax=(965.285,'K')), NASAPolynomial(coeffs=[14.1507,0.0381753,-1.33221e-05,2.35604e-09,-1.64446e-13,80112.9,-38.0542], Tmin=(965.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([CH][CH2])C1(9221)',
    structure = SMILES('[CH2][CH]C1CC([CH][CH2])C1'),
    E0 = (700.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900779,0.053868,1.31645e-05,-4.44949e-08,1.75721e-11,84343.8,35.3269], Tmin=(100,'K'), Tmax=(1070.77,'K')), NASAPolynomial(coeffs=[12.2689,0.0438946,-1.83828e-05,3.48666e-09,-2.4717e-13,80046.5,-28.9931], Tmin=(1070.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC[C]([CH2])C=C(5691)',
    structure = SMILES('[CH2]C=CCC([CH2])=C[CH2]'),
    E0 = (432.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,764.157],'cm^-1')),
        HinderedRotor(inertia=(0.0205827,'amu*angstrom^2'), symmetry=1, barrier=(21.2846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925992,'amu*angstrom^2'), symmetry=1, barrier=(21.2904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0205877,'amu*angstrom^2'), symmetry=1, barrier=(21.2888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.93926,'amu*angstrom^2'), symmetry=1, barrier=(90.5713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925935,'amu*angstrom^2'), symmetry=1, barrier=(21.2891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0502571,0.0722772,-2.38551e-05,-2.10775e-08,1.32161e-11,52129.6,32.212], Tmin=(100,'K'), Tmax=(1011.58,'K')), NASAPolynomial(coeffs=[17.8208,0.0334307,-1.28457e-05,2.37329e-09,-1.68106e-13,46926.6,-61.6719], Tmin=(1011.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])C=C(5695)',
    structure = SMILES('[CH2][C]=CCC([CH2])C=C'),
    E0 = (590.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0521473,0.0771334,-5.17676e-05,1.46563e-08,-4.82521e-13,71152.5,35.5549], Tmin=(100,'K'), Tmax=(1102.42,'K')), NASAPolynomial(coeffs=[15.0577,0.0358698,-1.3559e-05,2.39716e-09,-1.62227e-13,67043,-41.9342], Tmin=(1102.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CCC=C[CH2](5277)',
    structure = SMILES('[CH2]C=CCC=C[CH2]'),
    E0 = (319.697,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,654.584,654.666],'cm^-1')),
        HinderedRotor(inertia=(0.0715969,'amu*angstrom^2'), symmetry=1, barrier=(1.64714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02589,'amu*angstrom^2'), symmetry=1, barrier=(23.6022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02592,'amu*angstrom^2'), symmetry=1, barrier=(23.602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02593,'amu*angstrom^2'), symmetry=1, barrier=(23.6024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891956,0.0548009,-2.00953e-06,-3.44973e-08,1.64357e-11,38574.5,27.4605], Tmin=(100,'K'), Tmax=(1007.62,'K')), NASAPolynomial(coeffs=[14.652,0.0308083,-1.18922e-05,2.21117e-09,-1.57404e-13,34246.5,-46.7442], Tmin=(1007.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]CC=C(4889)',
    structure = SMILES('[CH2]C=CCC=C'),
    E0 = (204.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,272.006,272.35],'cm^-1')),
        HinderedRotor(inertia=(0.411509,'amu*angstrom^2'), symmetry=1, barrier=(21.5865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00227176,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409484,'amu*angstrom^2'), symmetry=1, barrier=(21.5978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56194,0.0414538,1.08814e-05,-4.06106e-08,1.75672e-11,24660.9,23.662], Tmin=(100,'K'), Tmax=(1001.11,'K')), NASAPolynomial(coeffs=[12.2191,0.0272829,-1.04542e-05,1.94474e-09,-1.38777e-13,21103.4,-34.8748], Tmin=(1001.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH][CH2](5531)',
    structure = SMILES('[CH2][CH][CH][CH2]'),
    E0 = (655.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1802.64],'cm^-1')),
        HinderedRotor(inertia=(0.00215831,'amu*angstrom^2'), symmetry=1, barrier=(4.96293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00215621,'amu*angstrom^2'), symmetry=1, barrier=(4.95965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214837,'amu*angstrom^2'), symmetry=1, barrier=(4.95028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82147,0.0344596,-5.62911e-05,6.39381e-08,-2.62695e-11,78855.1,21.3291], Tmin=(100,'K'), Tmax=(865.068,'K')), NASAPolynomial(coeffs=[-1.36886,0.0321046,-1.45274e-05,2.71435e-09,-1.84229e-13,80393.2,45.6373], Tmin=(865.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])CC=C[CH2](9222)',
    structure = SMILES('[CH2][CH][C]([CH2])CC=C[CH2]'),
    E0 = (812.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,232.449,628.887,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00259884,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00259884,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00259884,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00259884,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00259884,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00259884,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37302,0.0686926,-1.57463e-05,-8.87497e-08,9.19151e-11,97769.7,35.5106], Tmin=(100,'K'), Tmax=(475.815,'K')), NASAPolynomial(coeffs=[5.51715,0.0524984,-2.34696e-05,4.42158e-09,-3.06394e-13,97164.3,16.3778], Tmin=(475.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(812.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C=C(3739)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,230.296,230.449,230.657],'cm^-1')),
        HinderedRotor(inertia=(1.34138,'amu*angstrom^2'), symmetry=1, barrier=(50.5187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,44095,-3.44607], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH2](499)',
    structure = SMILES('[CH2][CH]C([CH2])[CH2]'),
    E0 = (636.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1200.87],'cm^-1')),
        HinderedRotor(inertia=(0.112344,'amu*angstrom^2'), symmetry=1, barrier=(2.58301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112403,'amu*angstrom^2'), symmetry=1, barrier=(2.58438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252635,'amu*angstrom^2'), symmetry=1, barrier=(2.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252511,'amu*angstrom^2'), symmetry=1, barrier=(2.5843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95305,0.0434597,-3.02543e-05,1.30344e-08,-2.49244e-12,76589.1,26.1277], Tmin=(100,'K'), Tmax=(1189.32,'K')), NASAPolynomial(coeffs=[6.71646,0.027439,-1.00484e-05,1.70807e-09,-1.11583e-13,75456,2.32115], Tmin=(1189.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH]C=C[CH2](9223)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]C=C[CH2]'),
    E0 = (767.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,241.65,831.818,1933.2],'cm^-1')),
        HinderedRotor(inertia=(0.0974324,'amu*angstrom^2'), symmetry=1, barrier=(3.27035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974324,'amu*angstrom^2'), symmetry=1, barrier=(3.27035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974324,'amu*angstrom^2'), symmetry=1, barrier=(3.27035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974324,'amu*angstrom^2'), symmetry=1, barrier=(3.27035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974324,'amu*angstrom^2'), symmetry=1, barrier=(3.27035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974324,'amu*angstrom^2'), symmetry=1, barrier=(3.27035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.152487,0.0752921,-5.13771e-05,1.82244e-08,-2.65215e-12,92503.2,37.4613], Tmin=(100,'K'), Tmax=(1587.44,'K')), NASAPolynomial(coeffs=[15.9984,0.0353642,-1.36489e-05,2.38005e-09,-1.56907e-13,87472.2,-46.3091], Tmin=(1587.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[C]=C[CH2](9224)',
    structure = SMILES('[CH2][CH]C([CH2])C[C]=C[CH2]'),
    E0 = (864.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3010,987.5,1337.5,450,1655,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494846,0.0785656,-6.32465e-05,2.9139e-08,-5.79131e-12,104116,37.8004], Tmin=(100,'K'), Tmax=(1155.87,'K')), NASAPolynomial(coeffs=[10.4052,0.0442694,-1.87391e-05,3.46839e-09,-2.39037e-13,101825,-11.447], Tmin=(1155.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(864.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])[CH][CH2](9225)',
    structure = SMILES('[CH2][C]=CCC([CH2])[CH][CH2]'),
    E0 = (864.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494846,0.0785656,-6.32465e-05,2.9139e-08,-5.79131e-12,104116,37.8004], Tmin=(100,'K'), Tmax=(1155.87,'K')), NASAPolynomial(coeffs=[10.4052,0.0442694,-1.87391e-05,3.46839e-09,-2.39037e-13,101825,-11.447], Tmin=(1155.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(864.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[C]([CH2])C[CH2](5668)',
    structure = SMILES('[CH2]C=CC[C]([CH2])C[CH2]'),
    E0 = (617.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,206.938,813.876,1637.15],'cm^-1')),
        HinderedRotor(inertia=(0.148221,'amu*angstrom^2'), symmetry=1, barrier=(3.56797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148221,'amu*angstrom^2'), symmetry=1, barrier=(3.56797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148221,'amu*angstrom^2'), symmetry=1, barrier=(3.56797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148221,'amu*angstrom^2'), symmetry=1, barrier=(3.56797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148221,'amu*angstrom^2'), symmetry=1, barrier=(3.56797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148221,'amu*angstrom^2'), symmetry=1, barrier=(3.56797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.473752,0.081835,-7.05277e-05,3.84045e-08,-9.39963e-12,74411.6,35.9936], Tmin=(100,'K'), Tmax=(939.868,'K')), NASAPolynomial(coeffs=[7.5847,0.0515715,-2.22284e-05,4.14505e-09,-2.86836e-13,73075,2.12834], Tmin=(939.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(RCCJ) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][C](C)CC=C[CH2](5798)',
    structure = SMILES('[CH2][CH][C](C)CC=C[CH2]'),
    E0 = (607.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,184.098,405.743,2978.29],'cm^-1')),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793825,0.077034,-6.20887e-05,3.37689e-08,-9.00501e-12,73130.4,35.6306], Tmin=(100,'K'), Tmax=(828.936,'K')), NASAPolynomial(coeffs=[4.87591,0.0573368,-2.64468e-05,5.10518e-09,-3.60577e-13,72453.6,16.7027], Tmin=(828.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C[CH][C]C(5864)',
    structure = SMILES('[CH2]C(C=C)C[CH][C]C'),
    E0 = (728.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.43916,0.0801256,-6.29041e-05,2.92675e-08,-5.98444e-12,87701.8,35.8475], Tmin=(100,'K'), Tmax=(1116.03,'K')), NASAPolynomial(coeffs=[9.33247,0.048251,-2.00632e-05,3.67644e-09,-2.51863e-13,85716.7,-8.03392], Tmin=(1116.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(728.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])C[CH2](5670)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])C[CH2]'),
    E0 = (573.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0533807,0.0756462,-3.80211e-05,-9.44003e-10,4.98608e-12,69110.7,35.2471], Tmin=(100,'K'), Tmax=(1054.59,'K')), NASAPolynomial(coeffs=[14.9698,0.0392785,-1.50385e-05,2.69876e-09,-1.85158e-13,64840.7,-42.8372], Tmin=(1054.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C=C[CH2](5799)',
    structure = SMILES('[CH2][CH]C(C)[CH]C=C[CH2]'),
    E0 = (562.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162606,0.0732669,-3.76535e-05,4.08843e-09,1.69225e-12,67838.6,35.6441], Tmin=(100,'K'), Tmax=(1226.2,'K')), NASAPolynomial(coeffs=[14.9314,0.0406437,-1.67731e-05,3.08096e-09,-2.11488e-13,63047.4,-43.3863], Tmin=(1226.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C]([CH2])[CH]C(5852)',
    structure = SMILES('[CH2]C=CC[C]([CH2])[CH]C'),
    E0 = (606.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,269.319,931.293,2084.77],'cm^-1')),
        HinderedRotor(inertia=(0.105449,'amu*angstrom^2'), symmetry=1, barrier=(3.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105449,'amu*angstrom^2'), symmetry=1, barrier=(3.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105449,'amu*angstrom^2'), symmetry=1, barrier=(3.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105449,'amu*angstrom^2'), symmetry=1, barrier=(3.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105449,'amu*angstrom^2'), symmetry=1, barrier=(3.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105449,'amu*angstrom^2'), symmetry=1, barrier=(3.3376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596437,0.0802031,-7.24225e-05,4.54783e-08,-1.3376e-11,73118.6,36.3259], Tmin=(100,'K'), Tmax=(782.502,'K')), NASAPolynomial(coeffs=[5.49147,0.0551812,-2.44586e-05,4.61557e-09,-3.21213e-13,72352.6,13.9107], Tmin=(782.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C[C][CH]C(5863)',
    structure = SMILES('[CH2]C(C=C)C[C][CH]C'),
    E0 = (728.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282666,0.0811219,-6.31384e-05,2.82724e-08,-5.43887e-12,87708,36.576], Tmin=(100,'K'), Tmax=(1197.14,'K')), NASAPolynomial(coeffs=[10.9935,0.0453337,-1.82964e-05,3.30068e-09,-2.24011e-13,85143.5,-17.025], Tmin=(1197.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(728.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])[CH]C(5672)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])[CH]C'),
    E0 = (562.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.138701,0.0743857,-4.08173e-05,6.53511e-09,1.25846e-12,67819.4,35.7176], Tmin=(100,'K'), Tmax=(1170.17,'K')), NASAPolynomial(coeffs=[14.1678,0.0407498,-1.60568e-05,2.8864e-09,-1.96273e-13,63555.7,-38.3585], Tmin=(1170.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])C[CH2](5675)',
    structure = SMILES('[CH2]C=[C]CC([CH2])C[CH2]'),
    E0 = (670.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3010,987.5,1337.5,450,1655,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0027103,0.083425,-6.46211e-05,2.7435e-08,-4.85664e-12,80740.7,37.0269], Tmin=(100,'K'), Tmax=(1317.05,'K')), NASAPolynomial(coeffs=[14.0282,0.0408107,-1.60862e-05,2.86695e-09,-1.93073e-13,77044.8,-34.5281], Tmin=(1317.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C(C)C[C]=C[CH2](5801)',
    structure = SMILES('[CH2][CH]C(C)C[C]=C[CH2]'),
    E0 = (659.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.458677,0.0771922,-5.20348e-05,1.82049e-08,-2.68271e-12,79452.5,36.1384], Tmin=(100,'K'), Tmax=(1508.7,'K')), NASAPolynomial(coeffs=[13.1465,0.0435532,-1.85898e-05,3.4262e-09,-2.33795e-13,75624.1,-30.2907], Tmin=(1508.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH][CH][CH]C)C=C(5873)',
    structure = SMILES('[CH2]C([CH][CH][CH]C)C=C'),
    E0 = (668.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,336.587,450.931,1997.38,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0146901,'amu*angstrom^2'), symmetry=1, barrier=(2.37737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0146901,'amu*angstrom^2'), symmetry=1, barrier=(2.37737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0146901,'amu*angstrom^2'), symmetry=1, barrier=(2.37737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0146901,'amu*angstrom^2'), symmetry=1, barrier=(2.37737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0146901,'amu*angstrom^2'), symmetry=1, barrier=(2.37737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0146901,'amu*angstrom^2'), symmetry=1, barrier=(2.37737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745128,0.0777165,-8.04968e-05,6.56848e-08,-2.40392e-11,80563.1,40.1336], Tmin=(100,'K'), Tmax=(803.77,'K')), NASAPolynomial(coeffs=[2.70691,0.0580793,-2.54225e-05,4.72068e-09,-3.23383e-13,80566.7,33.0818], Tmin=(803.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=[C]CC([CH2])[CH]C(5677)',
    structure = SMILES('[CH2]C=[C]CC([CH2])[CH]C'),
    E0 = (659.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355652,0.0790998,-5.74216e-05,2.29751e-08,-3.92869e-12,79437.2,36.5074], Tmin=(100,'K'), Tmax=(1325.77,'K')), NASAPolynomial(coeffs=[11.7881,0.0446072,-1.83965e-05,3.35141e-09,-2.2829e-13,76405.8,-21.8714], Tmin=(1325.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])C[CH2](5679)',
    structure = SMILES('[CH2][C]=CCC([CH2])C[CH2]'),
    E0 = (670.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0027103,0.083425,-6.46211e-05,2.7435e-08,-4.85664e-12,80740.7,37.0269], Tmin=(100,'K'), Tmax=(1317.05,'K')), NASAPolynomial(coeffs=[14.0282,0.0408107,-1.60862e-05,2.86695e-09,-1.93073e-13,77044.8,-34.5281], Tmin=(1317.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC(C)[CH][CH2](5802)',
    structure = SMILES('[CH2][C]=CCC(C)[CH][CH2]'),
    E0 = (659.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.458677,0.0771922,-5.20348e-05,1.82049e-08,-2.68271e-12,79452.5,36.1384], Tmin=(100,'K'), Tmax=(1508.7,'K')), NASAPolynomial(coeffs=[13.1465,0.0435532,-1.85898e-05,3.4262e-09,-2.33795e-13,75624.1,-30.2907], Tmin=(1508.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C](C=C)C[CH][CH]C(5872)',
    structure = SMILES('[CH2]C=C([CH2])C[CH][CH]C'),
    E0 = (552.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.575,829.175,1950.73],'cm^-1')),
        HinderedRotor(inertia=(0.09477,'amu*angstrom^2'), symmetry=1, barrier=(3.23892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09477,'amu*angstrom^2'), symmetry=1, barrier=(3.23892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09477,'amu*angstrom^2'), symmetry=1, barrier=(3.23892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09477,'amu*angstrom^2'), symmetry=1, barrier=(3.23892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09477,'amu*angstrom^2'), symmetry=1, barrier=(3.23892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09477,'amu*angstrom^2'), symmetry=1, barrier=(3.23892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811185,0.0751994,-5.30189e-05,2.25294e-08,-4.48413e-12,66589.8,34.4462], Tmin=(100,'K'), Tmax=(1082.93,'K')), NASAPolynomial(coeffs=[6.55367,0.0539884,-2.36387e-05,4.44246e-09,-3.08656e-13,65346,6.28461], Tmin=(1082.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CCC([CH2])[CH]C(5681)',
    structure = SMILES('[CH2][C]=CCC([CH2])[CH]C'),
    E0 = (659.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355652,0.0790998,-5.74216e-05,2.29751e-08,-3.92869e-12,79437.2,36.5074], Tmin=(100,'K'), Tmax=(1325.77,'K')), NASAPolynomial(coeffs=[11.7881,0.0446072,-1.83965e-05,3.35141e-09,-2.2829e-13,76405.8,-21.8714], Tmin=(1325.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])C=C(5773)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])C=C'),
    E0 = (687.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,258.129,921.441,2070.15],'cm^-1')),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3541.48,'J/mol'), sigma=(6.54527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.17 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105516,0.0799229,-6.00862e-05,2.53045e-08,-4.46691e-12,82891.1,40.4157], Tmin=(100,'K'), Tmax=(1320.6,'K')), NASAPolynomial(coeffs=[13.1623,0.0403744,-1.51647e-05,2.62703e-09,-1.73834e-13,79442.6,-26.2068], Tmin=(1320.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C2H3(60)',
    structure = SMILES('[CH]=C'),
    E0 = (286.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,728.586,728.586,3521.64],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8307,-0.00147628,3.09002e-05,-3.80476e-08,1.43171e-11,34524.2,5.61949], Tmin=(100,'K'), Tmax=(933.662,'K')), NASAPolynomial(coeffs=[5.36086,0.00527345,-1.3196e-06,2.21564e-10,-1.68768e-14,33658.6,-4.76326], Tmin=(933.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]CC([CH2])[CH][CH2](9036)',
    structure = SMILES('[CH]CC([CH2])[CH][CH2]'),
    E0 = (858.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,194.818,1144.12,1359.09,1430.84,2179.2],'cm^-1')),
        HinderedRotor(inertia=(0.109694,'amu*angstrom^2'), symmetry=1, barrier=(2.6314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109694,'amu*angstrom^2'), symmetry=1, barrier=(2.6314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109694,'amu*angstrom^2'), symmetry=1, barrier=(2.6314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109694,'amu*angstrom^2'), symmetry=1, barrier=(2.6314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109694,'amu*angstrom^2'), symmetry=1, barrier=(2.6314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38929,0.057802,-4.47644e-05,2.00355e-08,-3.88638e-12,103392,30.313], Tmin=(100,'K'), Tmax=(1179.59,'K')), NASAPolynomial(coeffs=[8.58265,0.0334093,-1.37461e-05,2.50503e-09,-1.71017e-13,101695,-5.57891], Tmin=(1179.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(858.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC(C=C)C1[CH2](5651)',
    structure = SMILES('[CH2]C1CC(C=C)C1[CH2]'),
    E0 = (430.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.879583,0.0467983,5.53842e-05,-1.06872e-07,4.53358e-11,51926.8,31.1959], Tmin=(100,'K'), Tmax=(939.891,'K')), NASAPolynomial(coeffs=[17.1929,0.0324678,-9.67463e-06,1.64333e-09,-1.17147e-13,46426.7,-59.4414], Tmin=(939.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CC(C=C)C1(5687)',
    structure = SMILES('[CH2][CH]C1CC(C=C)C1'),
    E0 = (429.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00566,0.0460798,4.63464e-05,-8.62861e-08,3.41368e-11,51759.2,32.4972], Tmin=(100,'K'), Tmax=(998.961,'K')), NASAPolynomial(coeffs=[14.7115,0.0391593,-1.52768e-05,2.8986e-09,-2.10008e-13,46627.9,-45.5888], Tmin=(998.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)CCC=C(5715)',
    structure = SMILES('[CH2]C=C([CH2])CCC=C'),
    E0 = (291.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0321492,0.0750619,-2.95796e-05,-1.29433e-08,9.6101e-12,35230.4,32.64], Tmin=(100,'K'), Tmax=(1041.48,'K')), NASAPolynomial(coeffs=[16.8692,0.037449,-1.47263e-05,2.71742e-09,-1.90668e-13,30229.3,-56.6943], Tmin=(1041.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C=CC=C(9226)',
    structure = SMILES('[CH2]CC([CH2])C=CC=C'),
    E0 = (390.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0100195,0.0756202,-3.20655e-05,-1.3657e-08,1.13781e-11,47067.3,34.9178], Tmin=(100,'K'), Tmax=(978.67,'K')), NASAPolynomial(coeffs=[16.8982,0.0345938,-1.22236e-05,2.14463e-09,-1.47604e-13,42413,-53.1606], Tmin=(978.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)C=C(5699)',
    structure = SMILES('[CH2]C=C[CH]C(C)C=C'),
    E0 = (288.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347148,0.0646815,-2.36641e-06,-3.94338e-08,1.85634e-11,34847.8,31.1331], Tmin=(100,'K'), Tmax=(1015.38,'K')), NASAPolynomial(coeffs=[15.9939,0.0380906,-1.48601e-05,2.76368e-09,-1.96114e-13,29863.6,-53.4891], Tmin=(1015.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1C[CH][CH]CC1[CH2](9227)',
    structure = SMILES('[CH2]C1C[CH][CH]CC1[CH2]'),
    E0 = (566.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31692,0.0490347,1.56829e-05,-3.61766e-08,1.24724e-11,68285.1,29.8056], Tmin=(100,'K'), Tmax=(1145.19,'K')), NASAPolynomial(coeffs=[6.64231,0.0542384,-2.23128e-05,4.09364e-09,-2.81267e-13,65504.5,-3.42345], Tmin=(1145.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclohexane) + radical(cyclohexane) + radical(cyclohexane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1C[CH][CH]CC1(9190)',
    structure = SMILES('[CH2][CH]C1C[CH][CH]CC1'),
    E0 = (565.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63492,0.0456589,1.78013e-05,-3.27495e-08,9.752e-12,68109.8,31.1379], Tmin=(100,'K'), Tmax=(1319.93,'K')), NASAPolynomial(coeffs=[7.03705,0.0565746,-2.56127e-05,4.83985e-09,-3.34012e-13,64306.7,-5.42792], Tmin=(1319.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclohexane) + radical(cyclohexane) + radical(cyclohexane) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1[CH]CC[CH][CH]C1(9188)',
    structure = SMILES('[CH2]C1[CH]CC[CH][CH]C1'),
    E0 = (545.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2292,0.0158268,0.000181871,-2.72478e-07,1.14074e-10,65698.5,35.6359], Tmin=(100,'K'), Tmax=(907.385,'K')), NASAPolynomial(coeffs=[30.8008,0.00319117,8.14869e-06,-1.86027e-09,1.19745e-13,55485.5,-130.861], Tmin=(907.385,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(cycloheptane) + radical(cycloheptane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])C=C(5692)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])C=C'),
    E0 = (493.616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.173,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375906,0.0662644,-1.46302e-05,-2.70679e-08,1.4933e-11,59510.7,32.81], Tmin=(100,'K'), Tmax=(992.681,'K')), NASAPolynomial(coeffs=[15.6006,0.0350196,-1.29048e-05,2.32182e-09,-1.61998e-13,55004.9,-47.9991], Tmin=(992.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])C=C(5669)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])C=C'),
    E0 = (679.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,221.668,765.878,1131.5,3560.55],'cm^-1')),
        HinderedRotor(inertia=(0.0542449,'amu*angstrom^2'), symmetry=1, barrier=(1.80903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542449,'amu*angstrom^2'), symmetry=1, barrier=(1.80903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542449,'amu*angstrom^2'), symmetry=1, barrier=(1.80903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542449,'amu*angstrom^2'), symmetry=1, barrier=(1.80903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542449,'amu*angstrom^2'), symmetry=1, barrier=(1.80903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0542449,'amu*angstrom^2'), symmetry=1, barrier=(1.80903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773824,0.0749556,-5.57158e-05,2.57541e-08,-5.51742e-12,81862.1,39.2613], Tmin=(100,'K'), Tmax=(1033.32,'K')), NASAPolynomial(coeffs=[6.73913,0.0518639,-2.21954e-05,4.12782e-09,-2.8523e-13,80629.3,10.2866], Tmin=(1033.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])CC[C]=C(9228)',
    structure = SMILES('[CH2][CH]C([CH2])CC[C]=C'),
    E0 = (728.723,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524875,0.0798681,-6.35254e-05,3.00164e-08,-6.29334e-12,87767.3,37.8557], Tmin=(100,'K'), Tmax=(1082.17,'K')), NASAPolynomial(coeffs=[8.77882,0.0493586,-2.12353e-05,3.96317e-09,-2.74478e-13,85980.9,-2.61644], Tmin=(1082.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(728.723,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]CC[C]([CH2])C=C(5671)',
    structure = SMILES('[CH2][CH]CCC([CH2])=C[CH2]'),
    E0 = (563.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19506,0.0803016,-5.68951e-05,2.14725e-08,-3.40252e-12,67916.1,35.8683], Tmin=(100,'K'), Tmax=(1435.44,'K')), NASAPolynomial(coeffs=[13.7291,0.0425875,-1.74847e-05,3.16893e-09,-2.14706e-13,64030.7,-34.3176], Tmin=(1435.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CCCC([CH2])[CH][CH2](9229)',
    structure = SMILES('[CH]=CCCC([CH2])[CH][CH2]'),
    E0 = (737.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.357073,0.0803993,-6.14286e-05,2.63463e-08,-4.85144e-12,88888.9,38.3785], Tmin=(100,'K'), Tmax=(1239.5,'K')), NASAPolynomial(coeffs=[11.2461,0.0452593,-1.89034e-05,3.47408e-09,-2.38245e-13,86189.5,-16.4927], Tmin=(1239.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(737.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]CC([CH2])C[CH2](9230)',
    structure = SMILES('[CH]C=CCC([CH2])C[CH2]'),
    E0 = (651.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.177625,0.0835573,-5.54035e-05,1.95476e-08,-2.88253e-12,78506.1,38.0178], Tmin=(100,'K'), Tmax=(1544.06,'K')), NASAPolynomial(coeffs=[15.1099,0.0439538,-1.69301e-05,2.93629e-09,-1.92972e-13,73785.1,-42.3769], Tmin=(1544.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=CCC(C)[CH][CH2](5785)',
    structure = SMILES('[CH]C=CCC(C)[CH][CH2]'),
    E0 = (640.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.182605,0.0783223,-4.55027e-05,1.28639e-08,-1.47865e-12,77222.8,37.5066], Tmin=(100,'K'), Tmax=(1926.76,'K')), NASAPolynomial(coeffs=[17.9152,0.0415092,-1.68437e-05,2.94786e-09,-1.92036e-13,70389.4,-59.6731], Tmin=(1926.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=CCC([CH2])[CH]C(5840)',
    structure = SMILES('[CH]C=CCC([CH2])[CH]C'),
    E0 = (640.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.187718,0.0791985,-4.82295e-05,1.52122e-08,-2.01762e-12,77202,37.4687], Tmin=(100,'K'), Tmax=(1659.17,'K')), NASAPolynomial(coeffs=[13.8557,0.0462472,-1.84393e-05,3.24222e-09,-2.14018e-13,72666.6,-35.3915], Tmin=(1659.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(640.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
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
    E0 = (626.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (626.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (786.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1088.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1083.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1189.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1092.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1081.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1057.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1081.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (632.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (634.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (634.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (649.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (649.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (651.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (651.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (651.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (979.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1011.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (973.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (852.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (688.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (701.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (695.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (701.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (700.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (651.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (813.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (741.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (762.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (626.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (930.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1012.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (979.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (1076.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1076.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (784.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (768.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (834.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (751.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (746.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (784.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (890.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (731.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (771.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (703.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (1087.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (801.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (714.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (699.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (715.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (734.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (845.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (1201.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (635.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (635.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (715.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (705.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (705.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (693.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (714.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (653.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (706.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (756.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (804.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (928.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (819.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (894.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (802.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (773.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (773.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C=C(3743)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C(C=C)C[CH]C=C(4934)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH][CH]CCC=C[CH2](5895)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][CH2](502)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]CC=C[CH2](5280)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH][CH]CC([CH2])C=C(5247)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C]C([CH2])CC=C[CH2](9211)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C([CH][CH2])CC=C[CH2](9212)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C=CCC([CH2])[CH][CH2](9213)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH][CH]C([CH2])CC=C[CH2](9214)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=CCC1CC1[CH2](5686)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C1CC=CCC1[CH2](5688)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['C=CC1C[CH][CH]CC1(5712)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_2H] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=CCC(=C)C[CH2](9215)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=CC[C](C)C=C(5697)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][C](C=C)CC=CC(5706)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]CC([CH2])CC=C=C(9014)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][C]=CCC(C)C=C(5705)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C[CH]C([CH2])[CH][CH2](9216)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CC[C]([CH2])[CH][CH2](9217)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][C]([CH2])C[CH]C[CH2](9218)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C([CH2])CC1[CH]C1(9219)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C1[CH]CC([CH2])C1[CH2](9206)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.5398e+06,'s^-1'), n=1.09629, Ea=(62.0258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C1C[CH]C([CH2])C1(9220)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.47079e+07,'s^-1'), n=0.909323, Ea=(74.2834,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C1[CH]CC([CH2])[CH]C1(5663)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.15136e+07,'s^-1'), n=0.706682, Ea=(68.453,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C1CC([CH2])C1[CH2](9209)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(74.7931,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 71.8 to 74.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C1CC([CH][CH2])C1(9221)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(73.4584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 70.0 to 73.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2]C=CC[C]([CH2])C=C(5691)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2][C]=CCC([CH2])C=C(5695)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(T)(20)', '[CH2]C=CCC=C[CH2](5277)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.84912e+14,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH][CH2](502)', 'C=C[CH]CC=C(4889)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.000715653,'m^3/(mol*s)'), n=2.74787, Ea=(77.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-CdHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from 72.5 to 77.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.9336e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])CC=C[CH2](9222)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C=C(3739)', '[CH2][CH]C([CH2])[CH2](499)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[CH]C=C[CH2](9223)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2e+07,'m^3/(mol*s)'), n=1.78837e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-3R!H-R_Ext-4R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-3R!H-R_Ext-4R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-3R!H-R_Ext-4R!H-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[C]=C[CH2](9224)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', '[CH2][C]=CCC([CH2])[CH][CH2](9225)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=CC[C]([CH2])C[CH2](5668)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH][C](C)CC=C[CH2](5798)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C(C=C)C[CH][C]C(5864)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=C[CH]C([CH2])C[CH2](5670)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C(C)[CH]C=C[CH2](5799)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 85 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=CC[C]([CH2])[CH]C(5852)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C=C)C[C][CH]C(5863)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=C[CH]C([CH2])[CH]C(5672)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(5.13003e+06,'s^-1'), n=1.69583, Ea=(104.251,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R4HJ_1;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C=[C]CC([CH2])C[CH2](5675)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.572e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][CH]C(C)C[C]=C[CH2](5801)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([CH][CH][CH]C)C=C(5873)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C=[C]CC([CH2])[CH]C(5677)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][C]=CCC([CH2])C[CH2](5679)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(760143,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;XH_out]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][C]=CCC(C)[CH][CH2](5802)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][C](C=C)C[CH][CH]C(5872)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SMSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2][C]=CCC([CH2])[CH]C(5681)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_4;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2][CH]C([CH2])C([CH2])C=C(5773)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C2H3(60)', '[CH]CC([CH2])[CH][CH2](9036)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C1CC(C=C)C1[CH2](5651)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C1CC(C=C)C1(5687)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][C](C=C)CCC=C(5715)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]CC([CH2])C=CC=C(9226)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C=C[CH]C(C)C=C(5699)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C1C[CH][CH]CC1[CH2](9227)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(5.52313e+07,'s^-1'), n=0.767814, Ea=(66.7931,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C1C[CH][CH]CC1(9190)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2]C1[CH]CC[CH][CH]C1(9188)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.2, Ea=(27.196,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 46 used for R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction65',
    reactants = ['H(3)', '[CH2]C=C[CH]C([CH2])C=C(5692)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.35e+08,'cm^3/(mol*s)'), n=1.64, Ea=(1.58992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2557 used for Cds-CsH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'butadiene13(1350)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(0.106847,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2][CH]C[CH]C([CH2])C=C(5669)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH2][CH]C([CH2])CC[C]=C(9228)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.9054e+11,'s^-1'), n=0.853, Ea=(200.196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]CC[C]([CH2])C=C(5671)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.32587e+10,'s^-1'), n=0.723333, Ea=(193.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH]=CCCC([CH2])[CH][CH2](9229)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(8.32e+10,'s^-1'), n=0.77, Ea=(268.194,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 195 used for R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH
Exact match found for rate rule [R3H_SD;C_rad_out_H/NonDeC;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]=C[CH]CC([CH2])C[CH2](9230)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]C=CCC(C)[CH][CH2](5785)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH]C=CCC([CH2])[CH]C(5840)'],
    products = ['[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2168',
    isomers = [
        '[CH2][CH]C([CH2])CC=C[CH2](5800)',
    ],
    reactants = [
        ('[CH2][CH]C=C(3743)', 'butadiene13(1350)'),
        ('[CH2][CH]C=C(3743)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2168',
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


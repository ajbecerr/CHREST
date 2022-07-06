species(
    label = '[CH2][CH][CH]CC[CH][CH]C(4586)',
    structure = SMILES('[CH2][CH][CH]CC[CH][CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,217.748,2872.67,2957.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674542,0.0897095,-0.00013367,1.39584e-07,-5.55415e-11,89541.8,42.901], Tmin=(100,'K'), Tmax=(857.285,'K')), NASAPolynomial(coeffs=[-4.6876,0.0737521,-3.40518e-05,6.36129e-09,-4.32043e-13,91967,76.7267], Tmin=(857.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = 'm1_allyl(186)',
    structure = SMILES('C=C[CH]C'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800698,'amu*angstrom^2'), symmetry=1, barrier=(1.84096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.062381,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C[CH]CC[CH][CH]C(4584)',
    structure = SMILES('[CH2]C=CCC[CH][CH]C'),
    E0 = (416.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,203.721,802.482,1211.16,1619.86],'cm^-1')),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3576.97,'J/mol'), sigma=(6.55299,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.71 K, Pc=28.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998225,0.0725476,-4.49469e-05,1.68524e-08,-3.19921e-12,50197.8,34.9872], Tmin=(100,'K'), Tmax=(1038.59,'K')), NASAPolynomial(coeffs=[4.31352,0.0597792,-2.65059e-05,5.0152e-09,-3.49866e-13,49509.2,18.8673], Tmin=(1038.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[CH]C(4576)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[CH]C'),
    E0 = (752.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,820.044,1553.04,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90564.5,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90974.1,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH][CH]C(4547)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH][CH]C'),
    E0 = (752.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,820.044,1553.04,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90564.5,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90974.1,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH]C(3874)',
    structure = SMILES('[CH][CH]C'),
    E0 = (522.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1852.85,1853.13,1853.16],'cm^-1')),
        HinderedRotor(inertia=(0.0369856,'amu*angstrom^2'), symmetry=1, barrier=(8.28142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370949,'amu*angstrom^2'), symmetry=1, barrier=(8.27684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19279,0.0163284,-1.76328e-06,-3.55309e-09,1.20622e-12,62877.3,14.3587], Tmin=(100,'K'), Tmax=(1426.65,'K')), NASAPolynomial(coeffs=[5.59458,0.0149965,-6.04303e-06,1.1011e-09,-7.44871e-14,61642.3,-0.00873183], Tmin=(1426.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][CH]C[CH2](509)',
    structure = SMILES('[CH2][CH][CH]C[CH2]'),
    E0 = (631.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1539.83],'cm^-1')),
        HinderedRotor(inertia=(0.00318581,'amu*angstrom^2'), symmetry=1, barrier=(5.35941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318763,'amu*angstrom^2'), symmetry=1, barrier=(5.35914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318741,'amu*angstrom^2'), symmetry=1, barrier=(5.35718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318288,'amu*angstrom^2'), symmetry=1, barrier=(5.35416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30707,0.0440519,-5.39188e-05,5.50378e-08,-2.25294e-11,75983.8,26.7965], Tmin=(100,'K'), Tmax=(841.968,'K')), NASAPolynomial(coeffs=[-0.503948,0.0407601,-1.83983e-05,3.43135e-09,-2.34064e-13,77047.1,43.3784], Tmin=(841.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH][CH2](5537)',
    structure = SMILES('[CH][CH][CH2]'),
    E0 = (727.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1609,1609.01,1609.02],'cm^-1')),
        HinderedRotor(inertia=(0.0337841,'amu*angstrom^2'), symmetry=1, barrier=(5.78124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00314687,'amu*angstrom^2'), symmetry=1, barrier=(5.78111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42424,0.0149081,-5.18362e-06,2.15834e-10,1.32851e-13,87550.9,15.303], Tmin=(100,'K'), Tmax=(1973.84,'K')), NASAPolynomial(coeffs=[8.10291,0.00913487,-3.61425e-06,6.37549e-10,-4.11102e-14,84981.5,-12.2801], Tmin=(1973.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C[CH][CH]C(123)',
    structure = SMILES('[CH2]C[CH][CH]C'),
    E0 = (426.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,2053.01],'cm^-1')),
        HinderedRotor(inertia=(0.00174415,'amu*angstrom^2'), symmetry=1, barrier=(5.21902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226977,'amu*angstrom^2'), symmetry=1, barrier=(5.21865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217106,'amu*angstrom^2'), symmetry=1, barrier=(64.9248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00174461,'amu*angstrom^2'), symmetry=1, barrier=(5.21934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.20251,0.0328813,-8.80397e-06,-8.44726e-10,4.254e-13,51258.1,21.7496], Tmin=(100,'K'), Tmax=(2245.78,'K')), NASAPolynomial(coeffs=[16.6245,0.0207134,-8.51707e-06,1.3975e-09,-8.32898e-14,42269.4,-60.4535], Tmin=(2245.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([532.913,1391.12,1391.12,2779.21,3448.45,3448.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.948,0.0008276,8.34932e-06,-9.82634e-09,3.80104e-12,16425.4,0.336655], Tmin=(100,'K'), Tmax=(660.467,'K')), NASAPolynomial(coeffs=[3.2217,0.00522646,-1.64125e-06,2.58225e-10,-1.62579e-14,16521.3,3.53938], Tmin=(660.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH][CH]CC[CH][CH][CH2](9242)',
    structure = SMILES('[CH][CH]CC[CH][CH][CH2]'),
    E0 = (1021.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180.354,798.464,838.401,1923.11,3012.99,3211.94,4000],'cm^-1')),
        HinderedRotor(inertia=(0.106846,'amu*angstrom^2'), symmetry=1, barrier=(2.45784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106846,'amu*angstrom^2'), symmetry=1, barrier=(2.45784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106846,'amu*angstrom^2'), symmetry=1, barrier=(2.45784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106846,'amu*angstrom^2'), symmetry=1, barrier=(2.45784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106846,'amu*angstrom^2'), symmetry=1, barrier=(2.45784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106846,'amu*angstrom^2'), symmetry=1, barrier=(2.45784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03987,0.0771788,-0.000114043,1.12873e-07,-4.32911e-11,122913,37.7863], Tmin=(100,'K'), Tmax=(857.161,'K')), NASAPolynomial(coeffs=[-0.440561,0.0552983,-2.53726e-05,4.7255e-09,-3.20435e-13,124224,50.8696], Tmin=(857.161,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1021.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH][CH]CC[C][CH]C(9420)',
    structure = SMILES('[CH2][CH][CH]CC[C][CH]C'),
    E0 = (997.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,219.919,1019.48,1181.09,1417.3,1784.09,1944.58],'cm^-1')),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356329,0.0924219,-0.000129661,1.23597e-07,-4.67871e-11,120077,41.3856], Tmin=(100,'K'), Tmax=(846.013,'K')), NASAPolynomial(coeffs=[0.570016,0.0635577,-2.90984e-05,5.43687e-09,-3.70331e-13,121038,46.2816], Tmin=(846.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(997.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH][C]C(9421)',
    structure = SMILES('[CH2][CH][CH]CC[CH][C]C'),
    E0 = (997.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,218.697,1141.3,1233.71,1428.42,1786.59,1932.4],'cm^-1')),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928716,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396181,0.0928662,-0.000134739,1.31829e-07,-5.0562e-11,120076,41.0707], Tmin=(100,'K'), Tmax=(848.024,'K')), NASAPolynomial(coeffs=[-0.517999,0.065501,-3.03036e-05,5.67994e-09,-3.87176e-13,121370,52.0461], Tmin=(848.024,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(997.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][C]CC[CH][CH]C(9422)',
    structure = SMILES('[CH2][CH][C]CC[CH][CH]C'),
    E0 = (997.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,219.919,1019.48,1181.09,1417.3,1784.09,1944.58],'cm^-1')),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0811966,'amu*angstrom^2'), symmetry=1, barrier=(2.60421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356329,0.0924219,-0.000129661,1.23597e-07,-4.67871e-11,120077,41.3856], Tmin=(100,'K'), Tmax=(846.013,'K')), NASAPolynomial(coeffs=[0.570016,0.0635577,-2.90984e-05,5.43687e-09,-3.70331e-13,121038,46.2816], Tmin=(846.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(997.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC[CH][CH]C(9423)',
    structure = SMILES('[CH2][C][CH]CC[CH][CH]C'),
    E0 = (997.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,218.697,1141.31,1233.71,1428.42,1786.59,1932.4],'cm^-1')),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0928718,'amu*angstrom^2'), symmetry=1, barrier=(2.92288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396121,0.092867,-0.000134742,1.31834e-07,-5.05642e-11,120076,41.0709], Tmin=(100,'K'), Tmax=(848.015,'K')), NASAPolynomial(coeffs=[-0.517874,0.0655008,-3.03034e-05,5.6799e-09,-3.87173e-13,121370,52.0454], Tmin=(848.015,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(997.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC[CH][CH]C(9424)',
    structure = SMILES('[CH][CH][CH]CC[CH][CH]C'),
    E0 = (986.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,186.018,436.859,470,2158.09,2565.71,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0799491,'amu*angstrom^2'), symmetry=1, barrier=(2.52427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.4314,0.0514131,-1.75357e-05,1.03806e-09,1.94735e-13,118581,28.2456], Tmin=(100,'K'), Tmax=(2792.21,'K')), NASAPolynomial(coeffs=[66.2181,-0.0117684,2.79768e-06,-5.67612e-10,4.75884e-14,74201.8,-350.967], Tmin=(2792.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(986.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CCC1[CH]C(5891)',
    structure = SMILES('[CH2][CH]C1CCC1[CH]C'),
    E0 = (495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865289,0.0531297,2.35527e-05,-5.65402e-08,2.18621e-11,59661.1,35.0531], Tmin=(100,'K'), Tmax=(1054.79,'K')), NASAPolynomial(coeffs=[12.3798,0.0463409,-1.92348e-05,3.64834e-09,-2.59318e-13,55180.6,-30.8362], Tmin=(1054.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C=CCC[CH][CH][CH]C(5908)',
    structure = SMILES('C=CCC[CH][CH][CH]C'),
    E0 = (471.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,208.962,666.525,815.213,3054,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0384586,'amu*angstrom^2'), symmetry=1, barrier=(1.11693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384586,'amu*angstrom^2'), symmetry=1, barrier=(1.11693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384586,'amu*angstrom^2'), symmetry=1, barrier=(1.11693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384586,'amu*angstrom^2'), symmetry=1, barrier=(1.11693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384586,'amu*angstrom^2'), symmetry=1, barrier=(1.11693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0384586,'amu*angstrom^2'), symmetry=1, barrier=(1.11693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81185,0.0586634,-2.30843e-05,2.98268e-09,-2.73343e-14,56744.3,31.0123], Tmin=(100,'K'), Tmax=(2619.53,'K')), NASAPolynomial(coeffs=[42.6604,0.0134167,-6.10888e-06,9.36106e-10,-4.90049e-14,30514.5,-209.825], Tmin=(2619.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C=C[CH]CC[CH]C(5900)',
    structure = SMILES('[CH2]C=C[CH]CC[CH]C'),
    E0 = (363.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0817584,0.0767014,-4.51066e-05,1.30975e-08,-1.54142e-12,43828.1,34.6245], Tmin=(100,'K'), Tmax=(1916.63,'K')), NASAPolynomial(coeffs=[18.6807,0.0378848,-1.47274e-05,2.5304e-09,-1.63066e-13,36698.8,-67.2045], Tmin=(1916.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[CH][CH][CH]C(5893)',
    structure = SMILES('[CH2]C=CC[CH][CH][CH]C'),
    E0 = (610.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,288.757,1278,1801.64,3619.33],'cm^-1')),
        HinderedRotor(inertia=(0.0842265,'amu*angstrom^2'), symmetry=1, barrier=(2.91154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842265,'amu*angstrom^2'), symmetry=1, barrier=(2.91154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842265,'amu*angstrom^2'), symmetry=1, barrier=(2.91154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842265,'amu*angstrom^2'), symmetry=1, barrier=(2.91154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842265,'amu*angstrom^2'), symmetry=1, barrier=(2.91154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842265,'amu*angstrom^2'), symmetry=1, barrier=(2.91154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73735,0.0578721,-2.44378e-05,3.94478e-09,-1.89486e-13,73497,31.1156], Tmin=(100,'K'), Tmax=(2765.97,'K')), NASAPolynomial(coeffs=[46.0237,0.00742123,-3.66583e-06,5.26021e-10,-2.39969e-14,44904.4,-230.155], Tmin=(2765.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C[CH][CH]C(5894)',
    structure = SMILES('[CH2]C=C[CH]C[CH][CH]C'),
    E0 = (557.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,245.799,857.211,1337.4,1844.1],'cm^-1')),
        HinderedRotor(inertia=(0.125739,'amu*angstrom^2'), symmetry=1, barrier=(3.49139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125739,'amu*angstrom^2'), symmetry=1, barrier=(3.49139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125739,'amu*angstrom^2'), symmetry=1, barrier=(3.49139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125739,'amu*angstrom^2'), symmetry=1, barrier=(3.49139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125739,'amu*angstrom^2'), symmetry=1, barrier=(3.49139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125739,'amu*angstrom^2'), symmetry=1, barrier=(3.49139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40166,0.0657818,-3.36033e-05,7.97104e-09,-7.60679e-13,67150.6,32.3576], Tmin=(100,'K'), Tmax=(2154.11,'K')), NASAPolynomial(coeffs=[13.2282,0.0438208,-1.83109e-05,3.23825e-09,-2.11404e-13,62055.5,-33.7742], Tmin=(2154.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(Allyl_S) + radical(RCCJC) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C(3856)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (450.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2810.32],'cm^-1')),
        HinderedRotor(inertia=(0.00169393,'amu*angstrom^2'), symmetry=1, barrier=(9.48541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010072,'amu*angstrom^2'), symmetry=1, barrier=(56.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100766,'amu*angstrom^2'), symmetry=1, barrier=(56.4312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79546,0.0335488,-4.51222e-05,5.0661e-08,-2.13657e-11,54172.4,20.3345], Tmin=(100,'K'), Tmax=(854.972,'K')), NASAPolynomial(coeffs=[-1.26875,0.0345785,-1.53757e-05,2.86235e-09,-1.94746e-13,55524.7,43.1492], Tmin=(854.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][CH][CH]C(9425)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH][CH]C'),
    E0 = (938.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3012.5,3025,3037.5,3050,390,398.75,407.5,416.25,425,1340,1345,1350,1355,1360,335,343.75,352.5,361.25,370,3000,3100,440,815,1455,1000,964.489,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932885,0.0901274,-0.000157371,1.75301e-07,-7.05875e-11,112915,44.5374], Tmin=(100,'K'), Tmax=(871.926,'K')), NASAPolynomial(coeffs=[-8.98133,0.0777523,-3.65489e-05,6.81925e-09,-4.60133e-13,116843,103.621], Tmin=(871.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(938.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C[CH][CH]C(9426)',
    structure = SMILES('[CH2][CH][CH][CH]C[CH][CH]C'),
    E0 = (938.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3012.5,3025,3037.5,3050,390,398.75,407.5,416.25,425,1340,1345,1350,1355,1360,335,343.75,352.5,361.25,370,3000,3100,440,815,1455,1000,964.489,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0098824,'amu*angstrom^2'), symmetry=1, barrier=(6.89931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932885,0.0901274,-0.000157371,1.75301e-07,-7.05875e-11,112915,44.5374], Tmin=(100,'K'), Tmax=(871.926,'K')), NASAPolynomial(coeffs=[-8.98133,0.0777523,-3.65489e-05,6.81925e-09,-4.60133e-13,116843,103.621], Tmin=(871.926,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(938.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH][CH][CH2](9026)',
    structure = SMILES('[CH2][CH][CH]CC[CH][CH][CH2]'),
    E0 = (948.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2672.96,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.036642,'amu*angstrom^2'), symmetry=1, barrier=(6.60198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702989,0.0905454,-0.000144462,1.52224e-07,-6.01034e-11,114225,43.8922], Tmin=(100,'K'), Tmax=(863.603,'K')), NASAPolynomial(coeffs=[-4.8421,0.0713671,-3.32299e-05,6.20584e-09,-4.20356e-13,116856,79.5166], Tmin=(863.603,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(948.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]C[CH]C(9427)',
    structure = SMILES('[CH2][CH][CH]C[CH]C[CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,187.672,2187.37,2892.19,3981.3,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674715,0.0897072,-0.000133661,1.39571e-07,-5.55352e-11,89541.8,42.9004], Tmin=(100,'K'), Tmax=(857.306,'K')), NASAPolynomial(coeffs=[-4.68798,0.0737528,-3.40522e-05,6.36139e-09,-4.32051e-13,91967.1,76.7288], Tmin=(857.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C[CH][CH]C(5887)',
    structure = SMILES('[CH2][CH]C[CH]C[CH][CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,217.748,2872.67,2957.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674542,0.0897095,-0.00013367,1.39584e-07,-5.55415e-11,89541.8,42.901], Tmin=(100,'K'), Tmax=(857.285,'K')), NASAPolynomial(coeffs=[-4.6876,0.0737521,-3.40518e-05,6.36129e-09,-4.32043e-13,91967,76.7267], Tmin=(857.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH]C[CH2](9096)',
    structure = SMILES('[CH2][CH][CH]CC[CH]C[CH2]'),
    E0 = (754.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,938.437,1411.63,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62777,0.0572323,-2.22985e-05,2.74783e-09,-1.63911e-14,90698.6,31.0442], Tmin=(100,'K'), Tmax=(2806,'K')), NASAPolynomial(coeffs=[60.9604,-0.00487899,4.17435e-07,-1.57664e-10,2.0493e-14,50800.5,-318.468], Tmin=(2806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]CC[CH]C(9428)',
    structure = SMILES('[CH2][CH][CH][CH]CC[CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,187.672,2187.37,2892.19,3981.3,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0304893,'amu*angstrom^2'), symmetry=1, barrier=(3.85225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674715,0.0897072,-0.000133661,1.39571e-07,-5.55352e-11,89541.8,42.9004], Tmin=(100,'K'), Tmax=(857.306,'K')), NASAPolynomial(coeffs=[-4.68798,0.0737528,-3.40522e-05,6.36139e-09,-4.32051e-13,91967.1,76.7288], Tmin=(857.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH][CH]C(5888)',
    structure = SMILES('[CH2][CH]CC[CH][CH][CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,217.748,2872.67,2957.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674542,0.0897095,-0.00013367,1.39584e-07,-5.55415e-11,89541.8,42.901], Tmin=(100,'K'), Tmax=(857.285,'K')), NASAPolynomial(coeffs=[-4.6876,0.0737521,-3.40518e-05,6.36129e-09,-4.32043e-13,91967,76.7267], Tmin=(857.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][CH]CC(9094)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH]CC'),
    E0 = (743.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,203.401,263.013,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716103,0.0901321,-0.000138665,1.47698e-07,-5.92618e-11,89540.5,42.5801], Tmin=(100,'K'), Tmax=(857.975,'K')), NASAPolynomial(coeffs=[-5.78083,0.0757047,-3.52624e-05,6.60569e-09,-4.49e-13,92301.2,82.5203], Tmin=(857.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C[CH][CH]C(9429)',
    structure = SMILES('[CH2]C[CH][CH]C[CH][CH]C'),
    E0 = (743.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,203.401,263.013,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716103,0.0901321,-0.000138665,1.47698e-07,-5.92618e-11,89540.5,42.5801], Tmin=(100,'K'), Tmax=(857.975,'K')), NASAPolynomial(coeffs=[-5.78083,0.0757047,-3.52624e-05,6.60569e-09,-4.49e-13,92301.2,82.5203], Tmin=(857.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CCC[CH][CH2](9430)',
    structure = SMILES('[CH2][CH][CH]CCC[CH][CH2]'),
    E0 = (754.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,213.636,754.186,2524.87,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413907,'amu*angstrom^2'), symmetry=1, barrier=(0.180735,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42237,0.0587917,-2.43146e-05,3.7173e-09,-1.59679e-13,90706.5,31.9467], Tmin=(100,'K'), Tmax=(2943.33,'K')), NASAPolynomial(coeffs=[66.4007,-0.0105801,2.77507e-06,-5.46706e-10,4.35076e-14,46609.2,-351.807], Tmin=(2943.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C[CH]CC(9095)',
    structure = SMILES('[CH2][CH][CH][CH]C[CH]CC'),
    E0 = (743.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,203.401,263.013,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716103,0.0901321,-0.000138665,1.47698e-07,-5.92618e-11,89540.5,42.5801], Tmin=(100,'K'), Tmax=(857.975,'K')), NASAPolynomial(coeffs=[-5.78083,0.0757047,-3.52624e-05,6.60569e-09,-4.49e-13,92301.2,82.5203], Tmin=(857.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[CH][CH][CH]C(5889)',
    structure = SMILES('[CH2]C[CH]C[CH][CH][CH]C'),
    E0 = (743.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,233.578,2668.45,3266.21,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0221334,'amu*angstrom^2'), symmetry=1, barrier=(3.44302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716069,0.0901326,-0.000138666,1.47701e-07,-5.9263e-11,89540.5,42.5802], Tmin=(100,'K'), Tmax=(857.972,'K')), NASAPolynomial(coeffs=[-5.78076,0.0757046,-3.52624e-05,6.60567e-09,-4.48999e-13,92301.2,82.5199], Tmin=(857.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C[CH][CH]C(4585)',
    structure = SMILES('C[CH][CH][CH]C[CH][CH]C'),
    E0 = (732.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3012.5,3025,3037.5,3050,390,398.75,407.5,416.25,425,1340,1345,1350,1355,1360,335,343.75,352.5,361.25,370,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901409,0.0893282,-0.00014671,1.62828e-07,-6.60925e-11,88231.8,42.8638], Tmin=(100,'K'), Tmax=(867.327,'K')), NASAPolynomial(coeffs=[-8.81198,0.0801113,-3.73554e-05,6.97101e-09,-4.71511e-13,91948.3,100.055], Tmin=(867.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJCC)"""),
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
    E0 = (743.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (743.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (909.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (909.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1209.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1209.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1213.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1209.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1209.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1209.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1209.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1198.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (751.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (807.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (807.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (830.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (777.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (841.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (813.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (761.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1105.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1149.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1149.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1160.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (880.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (880.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (895.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (859.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (859.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (890.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (917.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (895.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (929.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (896.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (819.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['[CH2][CH]C=C(3743)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['C=C[CH]CC[CH][CH]C(4584)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH]C(4576)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH]C(3874)', '[CH2][CH][CH]C[CH2](509)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C[CH][CH]C(123)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH][CH]CC[CH][CH][CH2](9242)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][CH][CH]CC[C][CH]C(9420)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][CH][CH]CC[CH][C]C(9421)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH][C]CC[CH][CH]C(9422)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C][CH]CC[CH][CH]C(9423)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH][CH][CH]CC[CH][CH]C(9424)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['[CH2][CH]C1CCC1[CH]C(5891)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['C=CCC[CH][CH][CH]C(5908)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['[CH2]C=C[CH]CC[CH]C(5900)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C=CC[CH][CH][CH]C(5893)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]C=C[CH]C[CH][CH]C(5894)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][CH][CH]CCC=C[CH2](5895)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'm1_allyl(186)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][CH][CH2](5531)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH][CH]C[CH][CH][CH]C(9425)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C[CH][CH]C(9426)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH][CH]CC[CH][CH][CH2](9026)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.53274e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH]C[CH]C[CH]C(9427)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C[CH]C[CH][CH]C(5887)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH][CH]CC[CH]C[CH2](9096)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH][CH]CC[CH]C(9428)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]CC[CH][CH][CH]C(5888)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]C[CH][CH]CC(9094)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C[CH][CH]C[CH][CH]C(9429)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]CCC[CH][CH2](9430)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH][CH][CH]C[CH]CC(9095)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.07601e+07,'s^-1'), n=1.21949, Ea=(185.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R4HJ_2;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C[CH]C[CH][CH][CH]C(5889)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(645397,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Y_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['C[CH][CH][CH]C[CH][CH]C(4585)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.32767e+06,'s^-1'), n=1.53625, Ea=(75.6258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2155',
    isomers = [
        '[CH2][CH][CH]CC[CH][CH]C(4586)',
    ],
    reactants = [
        ('[CH2][CH]C=C(3743)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2155',
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


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
    label = '[CH2]C(C=C)C[CH][CH]C(4572)',
    structure = SMILES('[CH2]C(C=C)C[CH][CH]C'),
    E0 = (474.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,400,1161.01,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3518.52,'J/mol'), sigma=(6.51243,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.58 K, Pc=28.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32824,0.067241,-1.17818e-05,-6.33616e-08,5.75461e-11,57141.7,35.6254], Tmin=(100,'K'), Tmax=(513.968,'K')), NASAPolynomial(coeffs=[4.21148,0.0582682,-2.48959e-05,4.62596e-09,-3.19747e-13,56667.5,21.9041], Tmin=(513.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3496.45,'J/mol'), sigma=(6.67923,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=546.14 K, Pc=26.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674542,0.0897095,-0.00013367,1.39584e-07,-5.55415e-11,89541.8,42.901], Tmin=(100,'K'), Tmax=(857.285,'K')), NASAPolynomial(coeffs=[-4.6876,0.0737521,-3.40518e-05,6.36129e-09,-4.32043e-13,91967,76.7267], Tmin=(857.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])[CH]C(4531)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[CH]C'),
    E0 = (760.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,347.706,3366.17,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418259,0.0816099,-6.41805e-05,3.1146e-08,-6.78753e-12,91578.7,41.3492], Tmin=(100,'K'), Tmax=(1042.54,'K')), NASAPolynomial(coeffs=[8.11405,0.0520824,-2.16959e-05,3.97833e-09,-2.72654e-13,89974.1,3.90099], Tmin=(1042.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH][CH]C(555)',
    structure = SMILES('[CH2][CH]C[CH][CH]C'),
    E0 = (596.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,3607.59,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79917,0.0586099,-7.80163e-05,8.14675e-08,-3.30993e-11,71839.7,31.8124], Tmin=(100,'K'), Tmax=(852.592,'K')), NASAPolynomial(coeffs=[-2.4828,0.0533525,-2.41736e-05,4.49758e-09,-3.05589e-13,73491.1,57.1902], Tmin=(852.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[CH][CH]C(4161)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH]C'),
    E0 = (767.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,3248.34,3721.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36377,0.0743961,-0.000117792,1.28512e-07,-5.18955e-11,92377.3,38.1829], Tmin=(100,'K'), Tmax=(865.662,'K')), NASAPolynomial(coeffs=[-5.72141,0.0655337,-3.03502e-05,5.65576e-09,-3.82639e-13,95162.7,80.3458], Tmin=(865.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH][CH]CC([CH2])[CH][CH2](9261)',
    structure = SMILES('[CH][CH]CC([CH2])[CH][CH2]'),
    E0 = (1029.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,185.694,884.919,991.706,1889.24,2444.97,4000],'cm^-1')),
        HinderedRotor(inertia=(0.101501,'amu*angstrom^2'), symmetry=1, barrier=(2.35315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101501,'amu*angstrom^2'), symmetry=1, barrier=(2.35315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101501,'amu*angstrom^2'), symmetry=1, barrier=(2.35315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101501,'amu*angstrom^2'), symmetry=1, barrier=(2.35315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101501,'amu*angstrom^2'), symmetry=1, barrier=(2.35315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101501,'amu*angstrom^2'), symmetry=1, barrier=(2.35315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848697,0.0738779,-8.19987e-05,6.24212e-08,-2.07174e-11,123934,37.2366], Tmin=(100,'K'), Tmax=(797.673,'K')), NASAPolynomial(coeffs=[5.85425,0.0446608,-1.93162e-05,3.56412e-09,-2.43329e-13,123267,15.0401], Tmin=(797.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1029.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][C]C([CH2])C[CH][CH]C(9440)',
    structure = SMILES('[CH2][C]C([CH2])C[CH][CH]C'),
    E0 = (1005.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,244.848,856.184,1141.58,1482.83,1839.72],'cm^-1')),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120082,'amu*angstrom^2'), symmetry=1, barrier=(3.44425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.186399,0.0913286,-0.000112018,9.34608e-08,-3.25498e-11,121085,40.7651], Tmin=(100,'K'), Tmax=(839.606,'K')), NASAPolynomial(coeffs=[4.80219,0.0557481,-2.4172e-05,4.43013e-09,-2.99731e-13,120789,22.1562], Tmin=(839.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[C][CH]C(9441)',
    structure = SMILES('[CH2][CH]C([CH2])C[C][CH]C'),
    E0 = (1005.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.833997,0.0790698,-4.87515e-05,-2.02743e-08,3.6921e-11,121070,38.5545], Tmin=(100,'K'), Tmax=(540.398,'K')), NASAPolynomial(coeffs=[6.46446,0.053687,-2.35224e-05,4.39588e-09,-3.03656e-13,120223,12.6548], Tmin=(540.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH][C]C(9442)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH][C]C'),
    E0 = (1005.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,244.952,856.318,1141.76,1483.17,1840.29],'cm^-1')),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118114,'amu*angstrom^2'), symmetry=1, barrier=(3.45947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34458,0.0709968,-1.91605e-06,-1.38965e-07,1.39435e-10,121049,36.701], Tmin=(100,'K'), Tmax=(454.818,'K')), NASAPolynomial(coeffs=[5.78455,0.0548582,-2.42478e-05,4.51938e-09,-3.10182e-13,120408,16.1738], Tmin=(454.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C[CH][CH]C(9443)',
    structure = SMILES('[CH]C([CH][CH2])C[CH][CH]C'),
    E0 = (995.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180.746,372.15,707.666,1911.72,2514.28,3528.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40976,0.087886,-0.000109436,9.81733e-08,-3.70162e-11,119808,41.8341], Tmin=(100,'K'), Tmax=(809.584,'K')), NASAPolynomial(coeffs=[2.46238,0.0608317,-2.79737e-05,5.28758e-09,-3.64684e-13,120030,35.7881], Tmin=(809.584,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(995.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C[CH][CH]C(9444)',
    structure = SMILES('[CH][CH]C([CH2])C[CH][CH]C'),
    E0 = (994.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,181.776,539.911,733.407,1747.3,2123.85,3570.52,4000],'cm^-1')),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106688,'amu*angstrom^2'), symmetry=1, barrier=(2.48495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344204,0.0893504,-0.000113018,9.97772e-08,-3.63546e-11,119790,42.0658], Tmin=(100,'K'), Tmax=(838.444,'K')), NASAPolynomial(coeffs=[2.88743,0.0590168,-2.61889e-05,4.84721e-09,-3.2948e-13,120004,34.0593], Tmin=(838.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(994.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC1C[CH][CH]C(5871)',
    structure = SMILES('[CH2]C1CC1C[CH][CH]C'),
    E0 = (499.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538544,0.0672025,-2.83655e-05,1.02794e-09,1.78971e-12,60156.8,36.3991], Tmin=(100,'K'), Tmax=(1234.06,'K')), NASAPolynomial(coeffs=[10.4709,0.046814,-1.7933e-05,3.14403e-09,-2.09396e-13,56806.4,-17.2498], Tmin=(1234.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC([CH]C)C1[CH2](5831)',
    structure = SMILES('[CH2]C1CC([CH]C)C1[CH2]'),
    E0 = (496.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.748867,0.0538163,3.22692e-05,-7.59884e-08,3.22177e-11,59828.2,33.7121], Tmin=(100,'K'), Tmax=(964.146,'K')), NASAPolynomial(coeffs=[14.3818,0.040428,-1.40668e-05,2.49312e-09,-1.74599e-13,55192.8,-41.9672], Tmin=(964.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CC([CH]C)C1(5849)',
    structure = SMILES('[CH2][CH]C1CC([CH]C)C1'),
    E0 = (495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865289,0.0531297,2.35527e-05,-5.65402e-08,2.18621e-11,59661.1,35.0531], Tmin=(100,'K'), Tmax=(1054.79,'K')), NASAPolynomial(coeffs=[12.3798,0.0463409,-1.92348e-05,3.64834e-09,-2.59318e-13,55180.6,-30.8362], Tmin=(1054.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(=C)C[CH][CH]C(9445)',
    structure = SMILES('[CH2]CC(=C)C[CH][CH]C'),
    E0 = (468.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.837,0.0594316,4.36239e-05,-2.14556e-07,1.95552e-10,56397,33.4154], Tmin=(100,'K'), Tmax=(416.142,'K')), NASAPolynomial(coeffs=[3.63124,0.060615,-2.70728e-05,5.12485e-09,-3.57086e-13,56088,24.4149], Tmin=(416.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.368,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(C)C[CH][CH]C(4542)',
    structure = SMILES('C=C[C](C)C[CH][CH]C'),
    E0 = (401.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,321.227,439.64,1270.67,3361.02],'cm^-1')),
        HinderedRotor(inertia=(0.0484521,'amu*angstrom^2'), symmetry=1, barrier=(2.43544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0484521,'amu*angstrom^2'), symmetry=1, barrier=(2.43544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0484521,'amu*angstrom^2'), symmetry=1, barrier=(2.43544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0484521,'amu*angstrom^2'), symmetry=1, barrier=(2.43544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0484521,'amu*angstrom^2'), symmetry=1, barrier=(2.43544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0484521,'amu*angstrom^2'), symmetry=1, barrier=(2.43544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40085,0.0655675,-3.08904e-05,6.57729e-09,-5.4374e-13,48345.4,31.6366], Tmin=(100,'K'), Tmax=(2682.37,'K')), NASAPolynomial(coeffs=[22.7516,0.0337288,-1.30859e-05,2.15221e-09,-1.31317e-13,36891.3,-92.435], Tmin=(2682.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][C](C=C)CC[CH]C(5878)',
    structure = SMILES('[CH2]C=C([CH2])CC[CH]C'),
    E0 = (358.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0741056,0.082138,-5.46863e-05,1.89135e-08,-2.70643e-12,43243.8,35.0588], Tmin=(100,'K'), Tmax=(1591.5,'K')), NASAPolynomial(coeffs=[16.0441,0.0416278,-1.65056e-05,2.92013e-09,-1.94142e-13,38113.3,-50.1921], Tmin=(1591.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C=C[CH]C(9446)',
    structure = SMILES('[CH2]CC([CH2])C=C[CH]C'),
    E0 = (420.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.10099,0.0743444,-3.14533e-05,-6.42971e-09,6.47034e-12,50750.6,35.4985], Tmin=(100,'K'), Tmax=(1058.82,'K')), NASAPolynomial(coeffs=[14.1468,0.042467,-1.6305e-05,2.92855e-09,-2.00835e-13,46588.7,-38.6748], Tmin=(1058.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C=C[CH]C(4543)',
    structure = SMILES('[CH2][CH]C(C)C=C[CH]C'),
    E0 = (410.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.227167,0.0717799,-3.05062e-05,-2.05415e-09,3.4179e-12,49477.8,35.8337], Tmin=(100,'K'), Tmax=(1206.4,'K')), NASAPolynomial(coeffs=[13.9265,0.0441183,-1.81956e-05,3.34607e-09,-2.30002e-13,44879.9,-38.1849], Tmin=(1206.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_S) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3599.97,'J/mol'), sigma=(6.5858,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.31 K, Pc=28.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.117574,0.0787344,-5.46928e-05,2.01205e-08,-3.07431e-12,75530.1,38.683], Tmin=(100,'K'), Tmax=(1501.78,'K')), NASAPolynomial(coeffs=[14.8195,0.0395759,-1.55808e-05,2.75806e-09,-1.84022e-13,71114.3,-38.2241], Tmin=(1501.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[CH][CH]C(4158)',
    structure = SMILES('[CH2]C=CC[CH][CH]C'),
    E0 = (440.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,295.384,982.151,2188.21],'cm^-1')),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90721,0.0555768,-2.65323e-05,5.65859e-09,-4.66342e-13,53019.7,29.4014], Tmin=(100,'K'), Tmax=(2690.43,'K')), NASAPolynomial(coeffs=[20.1883,0.0283967,-1.13781e-05,1.90342e-09,-1.17395e-13,43183.1,-76.8867], Tmin=(2690.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=CC(553)',
    structure = SMILES('[CH2][CH]CC=CC'),
    E0 = (323.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1782.86],'cm^-1')),
        HinderedRotor(inertia=(0.133732,'amu*angstrom^2'), symmetry=1, barrier=(3.07475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133557,'amu*angstrom^2'), symmetry=1, barrier=(3.07073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134216,'amu*angstrom^2'), symmetry=1, barrier=(3.0859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13443,'amu*angstrom^2'), symmetry=1, barrier=(3.09082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80558,0.0490136,-2.40298e-05,5.39389e-09,-4.72015e-13,38969.5,26.5298], Tmin=(100,'K'), Tmax=(2579.8,'K')), NASAPolynomial(coeffs=[18.9745,0.022393,-8.55148e-06,1.39402e-09,-8.43992e-14,30111.1,-72.5714], Tmin=(2579.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][C]([CH2])C[CH][CH]C(9447)',
    structure = SMILES('[CH2][CH][C]([CH2])C[CH][CH]C'),
    E0 = (937.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1987.56,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0729592,'amu*angstrom^2'), symmetry=1, barrier=(1.87307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.499648,0.0947272,-0.000152004,1.56192e-07,-6.03739e-11,112855,43.4227], Tmin=(100,'K'), Tmax=(868.078,'K')), NASAPolynomial(coeffs=[-3.36376,0.068846,-3.17997e-05,5.90833e-09,-3.98645e-13,115172,70.9949], Tmin=(868.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(937.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][CH][CH]C(9448)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH][CH]C'),
    E0 = (946.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2163.58,2984.83,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213878,'amu*angstrom^2'), symmetry=1, barrier=(8.02976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.623638,0.0858552,-0.000115703,1.11385e-07,-4.28911e-11,113953,44.3997], Tmin=(100,'K'), Tmax=(843.049,'K')), NASAPolynomial(coeffs=[-0.325052,0.0637533,-2.90438e-05,5.426e-09,-3.69992e-13,115059,54.4215], Tmin=(843.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(946.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[CH][CH2](9025)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[CH][CH2]'),
    E0 = (957.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2680.8,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465157,0.0878423,-0.000114739,1.05124e-07,-3.91147e-11,115248,44.2], Tmin=(100,'K'), Tmax=(844.606,'K')), NASAPolynomial(coeffs=[1.5906,0.0604829,-2.70259e-05,5.00868e-09,-3.40222e-13,115844,43.612], Tmin=(844.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(957.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])C[CH][CH]C(5868)',
    structure = SMILES('[CH2]C[C]([CH2])C[CH][CH]C'),
    E0 = (742.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2549.1,2711.53,3987.01,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260347,'amu*angstrom^2'), symmetry=1, barrier=(10.097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.309686,0.0959179,-0.000140406,1.37527e-07,-5.21532e-11,89467,41.5492], Tmin=(100,'K'), Tmax=(865.006,'K')), NASAPolynomial(coeffs=[-1.28822,0.067951,-3.05979e-05,5.64503e-09,-3.80442e-13,91066.1,56.6723], Tmin=(865.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(742.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C](C)C[CH][CH]C(4545)',
    structure = SMILES('[CH2][CH][C](C)C[CH][CH]C'),
    E0 = (732.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1830.52,2847.58,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08413,0.0545227,-1.89188e-05,1.05859e-09,2.41154e-13,88019.5,28.4524], Tmin=(100,'K'), Tmax=(2678.98,'K')), NASAPolynomial(coeffs=[58.2003,-0.00143803,-1.49401e-06,1.83565e-10,-1.80944e-16,50110.4,-302.591], Tmin=(2678.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH]C[CH]C(9449)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]C[CH]C'),
    E0 = (752.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1124,2369.39,3999.83,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597476,'amu*angstrom^2'), symmetry=1, barrier=(6.01757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451066,0.084354,-8.78942e-05,6.99508e-08,-2.5277e-11,90577,42.4604], Tmin=(100,'K'), Tmax=(776.311,'K')), NASAPolynomial(coeffs=[3.65373,0.0603149,-2.68819e-05,5.04924e-09,-3.48763e-13,90306.8,29.2831], Tmin=(776.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH]C[CH2](5792)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH]C[CH2]'),
    E0 = (762.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1007.55,2415.44,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331576,0.0867925,-9.20103e-05,7.18669e-08,-2.5212e-11,91870.3,41.9491], Tmin=(100,'K'), Tmax=(780.083,'K')), NASAPolynomial(coeffs=[4.49909,0.0589563,-2.60505e-05,4.87045e-09,-3.35455e-13,91416.9,24.1394], Tmin=(780.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[CH][CH][CH]C(5869)',
    structure = SMILES('[CH2]CC([CH2])[CH][CH][CH]C'),
    E0 = (752.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,820.041,1553.04,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286781,'amu*angstrom^2'), symmetry=1, barrier=(2.19905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449456,0.0868469,-0.00010335,9.16627e-08,-3.41882e-11,90564.5,42.4705], Tmin=(100,'K'), Tmax=(830.581,'K')), NASAPolynomial(coeffs=[1.69728,0.0629526,-2.78981e-05,5.17628e-09,-3.52934e-13,90974.1,40.3958], Tmin=(830.581,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][CH][CH]C(4546)',
    structure = SMILES('[CH2][CH]C(C)[CH][CH][CH]C'),
    E0 = (741.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,2156.61,2460.57,3867.98,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88834,0.0594493,-2.43178e-05,3.32932e-09,-5.99579e-14,89181.6,34.3076], Tmin=(100,'K'), Tmax=(2631.57,'K')), NASAPolynomial(coeffs=[45.7947,0.0102635,-5.42006e-06,8.56886e-10,-4.5003e-14,61048.2,-224.753], Tmin=(2631.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])CC[CH]C(9450)',
    structure = SMILES('[CH2][CH][C]([CH2])CC[CH]C'),
    E0 = (742.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.252845,0.0941655,-0.000127766,1.1974e-07,-4.50015e-11,89481.9,41.7454], Tmin=(100,'K'), Tmax=(845.698,'K')), NASAPolynomial(coeffs=[0.886995,0.0649216,-2.93475e-05,5.46119e-09,-3.71465e-13,90313.2,44.341], Tmin=(845.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(742.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][CH]CC(9121)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH]CC'),
    E0 = (752.099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,712.696,1296.47,3907.28,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0361088,'amu*angstrom^2'), symmetry=1, barrier=(1.73931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465094,0.0851339,-9.42991e-05,8.01355e-08,-2.99924e-11,90576.8,42.2362], Tmin=(100,'K'), Tmax=(793.284,'K')), NASAPolynomial(coeffs=[2.63709,0.0621291,-2.80092e-05,5.27333e-09,-3.63996e-13,90611.4,34.6509], Tmin=(793.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C[CH][CH]C(4574)',
    structure = SMILES('[CH2][C]([CH]C)C[CH][CH]C'),
    E0 = (732.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24434,0.0541381,-1.90351e-05,1.53347e-09,1.31205e-13,87988.9,27.8064], Tmin=(100,'K'), Tmax=(2822.37,'K')), NASAPolynomial(coeffs=[67.9679,-0.0118015,3.05647e-06,-6.25078e-10,5.13848e-14,42311.4,-362.938], Tmin=(2822.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])[CH][CH2](5791)',
    structure = SMILES('[CH2][CH]CCC([CH2])[CH][CH2]'),
    E0 = (762.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,208.296,1250.46,2018.91,4000],'cm^-1')),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10531,0.0737795,-2.35e-05,-6.32665e-08,6.5968e-11,91836.8,39.5049], Tmin=(100,'K'), Tmax=(503.992,'K')), NASAPolynomial(coeffs=[5.23185,0.0577035,-2.5282e-05,4.73755e-09,-3.28201e-13,91209.1,20.3232], Tmin=(503.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][CH][CH]C(4575)',
    structure = SMILES('[CH2]C([CH]C)[CH][CH][CH]C'),
    E0 = (741.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,1904.01,2561.33,3951.84,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04732,0.0590742,-2.44524e-05,3.81605e-09,-1.72319e-13,89151.2,33.6666], Tmin=(100,'K'), Tmax=(2868.23,'K')), NASAPolynomial(coeffs=[56.0077,-0.000598494,-6.63779e-07,1.10938e-11,9.03774e-15,52935.7,-287.812], Tmin=(2868.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C[CH]CC(9120)',
    structure = SMILES('[CH2][CH][C]([CH2])C[CH]CC'),
    E0 = (742.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2454.3,2949.76,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0246122,'amu*angstrom^2'), symmetry=1, barrier=(11.6832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292435,0.0946132,-0.000132858,1.27992e-07,-4.8786e-11,89480.7,41.4314], Tmin=(100,'K'), Tmax=(847.812,'K')), NASAPolynomial(coeffs=[-0.200423,0.0668639,-3.0552e-05,5.7041e-09,-3.88298e-13,90645.1,50.1022], Tmin=(847.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(742.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])C[CH2](9085)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])C[CH2]'),
    E0 = (762.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,187.146,1307.08,3863.53,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294331,0.0887901,-0.00010221,8.5135e-08,-3.02778e-11,91858.9,42.259], Tmin=(100,'K'), Tmax=(829.918,'K')), NASAPolynomial(coeffs=[3.60694,0.0596931,-2.58869e-05,4.76055e-09,-3.23299e-13,91761.3,29.6196], Tmin=(829.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)[CH][CH2](4548)',
    structure = SMILES('[CH2][CH][CH]CC(C)[CH][CH2]'),
    E0 = (752.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,880.259,1174.39,3820.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.707,0.0622921,-2.82776e-05,5.28345e-09,-3.58346e-13,90474.2,34.1373], Tmin=(100,'K'), Tmax=(3091.93,'K')), NASAPolynomial(coeffs=[61.417,-0.0058351,9.76926e-07,-2.05723e-10,1.93015e-14,50428.2,-321.425], Tmin=(3091.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    E0 = (752.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (752.003,'kJ/mol'),
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
    E0 = (917.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1210.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1205.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1214.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1221.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1217.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1217.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1217.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1206.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1206.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (757.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (760.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (760.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (774.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (774.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (815.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (815.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (815.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (772.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (888.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (846.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (861.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (881.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (761.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (813.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1105.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1149.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (1158.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1169.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (910.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (893.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (889.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (904.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (868.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (869.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (868.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (898.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (909.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (904.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (883.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (870.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (836.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (804.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (815.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH]C=C(3743)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C(C=C)C[CH][CH]C(4572)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH][CH]CC[CH][CH]C(4586)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C[CH][CH]C(555)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C[CH][CH]C(4161)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH]C(3874)', '[CH2][CH]C([CH2])[CH2](499)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH][CH]CC([CH2])[CH][CH2](9261)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]C([CH2])C[CH][CH]C(9440)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[C][CH]C(9441)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[CH][C]C(9442)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C([CH][CH2])C[CH][CH]C(9443)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH][CH]C([CH2])C[CH][CH]C(9444)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C1CC1C[CH][CH]C(5871)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C1CC([CH]C)C1[CH2](5831)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH]C1CC([CH]C)C1(5849)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]CC(=C)C[CH][CH]C(9445)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C=C(C)C[CH][CH]C(4542)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][C](C=C)CC[CH]C(5878)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]CC([CH2])C=C[CH]C(9446)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH]C(C)C=C[CH]C(4543)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C](C=C)C[CH][CH]C(5872)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C([CH][CH][CH]C)C=C(5873)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH]C([CH2])CC=C[CH2](5800)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(T)(20)', '[CH2]C=CC[CH][CH]C(4158)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH][CH2](502)', '[CH2][CH]CC=CC(553)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'm1_allyl(186)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH][CH2](5531)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C[CH][CH]C(9447)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[CH][CH][CH]C(9448)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[CH][CH2](9025)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C[C]([CH2])C[CH][CH]C(5868)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH][C](C)C[CH][CH]C(4545)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]C([CH2])[CH]C[CH]C(9449)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]C([CH2])C[CH]C[CH2](5792)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]CC([CH2])[CH][CH][CH]C(5869)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH]C(C)[CH][CH][CH]C(4546)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH][C]([CH2])CC[CH]C(9450)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])[CH][CH]CC(9121)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][C]([CH]C)C[CH][CH]C(4574)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]CCC([CH2])[CH][CH2](5791)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C([CH]C)[CH][CH][CH]C(4575)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH][C]([CH2])C[CH]CC(9120)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH][CH]CC([CH2])C[CH2](9085)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2778.79,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH]C(4576)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2149',
    isomers = [
        '[CH2][CH]C([CH2])C[CH][CH]C(4547)',
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
    label = 'PDepNetwork #2149',
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


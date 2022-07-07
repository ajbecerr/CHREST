species(
    label = '[CH2]C([O])C(=C)[CH]O[O](20472)',
    structure = SMILES('[CH2]C([O])C(=C)[CH]O[O]'),
    E0 = (401.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,916.311],'cm^-1')),
        HinderedRotor(inertia=(0.0340544,'amu*angstrom^2'), symmetry=1, barrier=(3.07764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.929837,'amu*angstrom^2'), symmetry=1, barrier=(21.3788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359981,'amu*angstrom^2'), symmetry=1, barrier=(21.3844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359349,'amu*angstrom^2'), symmetry=1, barrier=(21.3742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0859817,0.0760865,-7.36512e-05,3.52531e-08,-6.5779e-12,48406.8,34.4714], Tmin=(100,'K'), Tmax=(1310.73,'K')), NASAPolynomial(coeffs=[19.2616,0.0175676,-6.68204e-06,1.19097e-09,-8.1103e-14,43380,-63.2283], Tmin=(1310.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(CJCO)"""),
)

species(
    label = 'vinoxy(1351)',
    structure = SMILES('[CH2]C=O'),
    E0 = (6.53237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1127.76,1127.76,1127.77,1945],'cm^-1')),
        HinderedRotor(inertia=(0.00342653,'amu*angstrom^2'), symmetry=1, barrier=(3.09254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51396,0.0045106,2.69938e-05,-3.71579e-08,1.42465e-11,808.756,8.86322], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[6.49564,0.00770518,-2.52913e-06,4.69083e-10,-3.5128e-14,-479.661,-9.13857], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.53237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C=CO[O](16806)',
    structure = SMILES('C=C=CO[O]'),
    E0 = (250.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.895395,'amu*angstrom^2'), symmetry=1, barrier=(20.5869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2049,0.0363638,-3.70353e-05,1.95556e-08,-4.06022e-12,30187.6,17.323], Tmin=(100,'K'), Tmax=(1179.28,'K')), NASAPolynomial(coeffs=[9.983,0.00998085,-3.47674e-06,5.84074e-10,-3.83218e-14,28353.2,-21.4843], Tmin=(1179.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'C=C([CH]O[O])C[CH][O](20382)',
    structure = SMILES('C=C([CH]O[O])C[CH][O]'),
    E0 = (378.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,263.452,263.536,264.358,1655.34],'cm^-1')),
        HinderedRotor(inertia=(0.00688885,'amu*angstrom^2'), symmetry=1, barrier=(45.1789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917495,'amu*angstrom^2'), symmetry=1, barrier=(45.1777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275991,'amu*angstrom^2'), symmetry=1, barrier=(13.6879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4356.36,'J/mol'), sigma=(7.16244,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.45 K, Pc=26.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799187,0.0776825,-9.30708e-05,5.78414e-08,-1.12702e-11,45603.6,31.3517], Tmin=(100,'K'), Tmax=(620.58,'K')), NASAPolynomial(coeffs=[8.92717,0.0357517,-1.69998e-05,3.27829e-09,-2.29702e-13,44393.4,-5.60588], Tmin=(620.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])C([C]=C)O[O](21020)',
    structure = SMILES('[CH2]C([O])C([C]=C)O[O]'),
    E0 = (517.722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,233.558,233.56,233.561],'cm^-1')),
        HinderedRotor(inertia=(0.00309033,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00309028,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378778,'amu*angstrom^2'), symmetry=1, barrier=(14.6624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378779,'amu*angstrom^2'), symmetry=1, barrier=(14.6624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4373.92,'J/mol'), sigma=(7.1993,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=683.20 K, Pc=26.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375961,0.0834533,-0.000104117,6.92465e-08,-1.85007e-11,62394.9,33.4043], Tmin=(100,'K'), Tmax=(911.005,'K')), NASAPolynomial(coeffs=[12.9999,0.0280242,-1.28501e-05,2.45754e-09,-1.72183e-13,60094.8,-26.322], Tmin=(911.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cds_S)"""),
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
    label = '[CH2]C([O])C(=C)[CH][O](22887)',
    structure = SMILES('[CH2]C(=C[O])C([CH2])[O]'),
    E0 = (297.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,443.121,445.108,446.3],'cm^-1')),
        HinderedRotor(inertia=(0.164279,'amu*angstrom^2'), symmetry=1, barrier=(23.2132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165695,'amu*angstrom^2'), symmetry=1, barrier=(23.2126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166191,'amu*angstrom^2'), symmetry=1, barrier=(23.2128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509454,0.0621697,-2.83854e-05,-2.18119e-08,1.66915e-11,35948.4,28.6219], Tmin=(100,'K'), Tmax=(942.138,'K')), NASAPolynomial(coeffs=[21.5298,0.00921533,-1.85529e-06,3.00633e-10,-2.53249e-14,30376.9,-80.084], Tmin=(942.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36761e-06,-4.93131e-09,1.45956e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.7,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00338e-07,1.59031e-10,-1.14892e-14,-1048.44,6.08305], Tmin=(1087.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH]=C([CH2])C([CH2])[O](15788)',
    structure = SMILES('[CH]C(=C)C([CH2])[O]'),
    E0 = (584.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,453.915,453.915,453.915,453.915,453.915],'cm^-1')),
        HinderedRotor(inertia=(0.346495,'amu*angstrom^2'), symmetry=1, barrier=(50.6608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346495,'amu*angstrom^2'), symmetry=1, barrier=(50.6608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346494,'amu*angstrom^2'), symmetry=1, barrier=(50.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.75,'J/mol'), sigma=(6.3037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.14 K, Pc=33.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07394,0.0564806,-3.76085e-05,8.30249e-09,8.9676e-13,70381.2,26.1567], Tmin=(100,'K'), Tmax=(1072.46,'K')), NASAPolynomial(coeffs=[13.2472,0.0236782,-9.35364e-06,1.69429e-09,-1.16695e-13,67045.5,-36.8022], Tmin=(1072.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C(=C)[CH]O[O](21398)',
    structure = SMILES('[CH2]C=C([CH2])[CH]O[O]'),
    E0 = (421.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,896.755],'cm^-1')),
        HinderedRotor(inertia=(1.04239,'amu*angstrom^2'), symmetry=1, barrier=(23.9665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04234,'amu*angstrom^2'), symmetry=1, barrier=(23.9655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159378,'amu*angstrom^2'), symmetry=1, barrier=(90.9568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159389,'amu*angstrom^2'), symmetry=1, barrier=(90.9509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938748,0.0571864,-3.19379e-05,-3.56891e-09,6.44623e-12,50839.2,27.1506], Tmin=(100,'K'), Tmax=(1001.96,'K')), NASAPolynomial(coeffs=[15.7724,0.0193533,-7.31489e-06,1.34988e-09,-9.61817e-14,46793.2,-49.799], Tmin=(1001.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = 'C=C([CH][O])[CH]O[O](22888)',
    structure = SMILES('[CH2]C([CH]O[O])=C[O]'),
    E0 = (238.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,400.309,400.625],'cm^-1')),
        HinderedRotor(inertia=(0.284663,'amu*angstrom^2'), symmetry=1, barrier=(32.4526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285185,'amu*angstrom^2'), symmetry=1, barrier=(32.4459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284925,'amu*angstrom^2'), symmetry=1, barrier=(32.4445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0002,0.0524333,-1.94965e-05,-2.53944e-08,1.70984e-11,28854.3,26.1125], Tmin=(100,'K'), Tmax=(938.089,'K')), NASAPolynomial(coeffs=[19.6705,0.00685316,-1.02764e-06,1.50224e-10,-1.46489e-14,23854.1,-70.7485], Tmin=(938.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]O[O](21387)',
    structure = SMILES('[CH]O[O]'),
    E0 = (471.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.569664,'amu*angstrom^2'), symmetry=1, barrier=(13.0977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33696,0.0153854,-2.48163e-05,1.94592e-08,-5.80272e-12,56680.3,11.4536], Tmin=(100,'K'), Tmax=(834.588,'K')), NASAPolynomial(coeffs=[6.14825,0.00191155,-5.99818e-07,1.15174e-10,-8.22083e-15,56211.1,-1.6009], Tmin=(834.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C([O])[C]=C(2851)',
    structure = SMILES('[CH2]C([O])[C]=C'),
    E0 = (490.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,352.978,353.215],'cm^-1')),
        HinderedRotor(inertia=(0.00135594,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179525,'amu*angstrom^2'), symmetry=1, barrier=(15.8404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87237,0.0412989,-3.27604e-05,1.18345e-08,-1.29105e-12,59070.3,21.8995], Tmin=(100,'K'), Tmax=(1121.63,'K')), NASAPolynomial(coeffs=[11.3451,0.0142396,-5.56323e-06,1.01285e-09,-7.00338e-14,56522.5,-26.7732], Tmin=(1121.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S)"""),
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
    label = '[CH2]C(=[C]O[O])C([CH2])[O](22522)',
    structure = SMILES('[CH2]C(=[C]O[O])C([CH2])[O]'),
    E0 = (672.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,271.6,840.951],'cm^-1')),
        HinderedRotor(inertia=(0.0229541,'amu*angstrom^2'), symmetry=1, barrier=(11.516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500898,'amu*angstrom^2'), symmetry=1, barrier=(11.5166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00722642,'amu*angstrom^2'), symmetry=1, barrier=(11.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500867,'amu*angstrom^2'), symmetry=1, barrier=(11.5159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484422,0.0801008,-0.000103151,6.96189e-08,-1.866e-11,81009.3,35.2413], Tmin=(100,'K'), Tmax=(912.949,'K')), NASAPolynomial(coeffs=[13.3545,0.023714,-1.05097e-05,1.97209e-09,-1.36523e-13,78659.2,-25.678], Tmin=(912.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]C([O])C(=C)[CH]O[O](21139)',
    structure = SMILES('[CH]C([O])C(=C)[CH]O[O]'),
    E0 = (637.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.193259,0.0738673,-7.31295e-05,3.54189e-08,-6.65878e-12,76862.4,34.4566], Tmin=(100,'K'), Tmax=(1304.16,'K')), NASAPolynomial(coeffs=[19.2462,0.0154301,-5.91742e-06,1.06128e-09,-7.26577e-14,71892.7,-62.5227], Tmin=(1304.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])C(=C)C1OO1(22889)',
    structure = SMILES('[CH2]C([O])C(=C)C1OO1'),
    E0 = (247.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440414,0.0647164,-3.05939e-05,-1.94035e-08,1.58934e-11,29934.4,30.7314], Tmin=(100,'K'), Tmax=(932.764,'K')), NASAPolynomial(coeffs=[20.4888,0.0131614,-3.03724e-06,4.61118e-10,-3.36618e-14,24696.9,-72.6214], Tmin=(932.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(dioxirane) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'C=C([CH]O[O])C1CO1(22890)',
    structure = SMILES('C=C([CH]O[O])C1CO1'),
    E0 = (148.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.864446,0.0563647,-1.03495e-05,-3.99463e-08,2.44313e-11,18027.3,28.2745], Tmin=(100,'K'), Tmax=(883.602,'K')), NASAPolynomial(coeffs=[17.9312,0.0152739,-1.99466e-06,7.63166e-11,3.95858e-16,13599.3,-59.9407], Tmin=(883.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C1OC(O[O])C1=C(22791)',
    structure = SMILES('[CH2]C1OC(O[O])C1=C'),
    E0 = (194.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431569,0.0682387,-5.76893e-05,2.40003e-08,-3.93518e-12,23510.8,28.2181], Tmin=(100,'K'), Tmax=(1468.17,'K')), NASAPolynomial(coeffs=[18.2405,0.0197183,-8.11648e-06,1.4901e-09,-1.02102e-13,18281.6,-64.5384], Tmin=(1468.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C1C([O])CC1O[O](22719)',
    structure = SMILES('C=C1C([O])CC1O[O]'),
    E0 = (224.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998817,0.0540942,-1.82466e-05,-1.54507e-08,9.45462e-12,27091.9,27.5173], Tmin=(100,'K'), Tmax=(1044.79,'K')), NASAPolynomial(coeffs=[15.5032,0.0233866,-9.79707e-06,1.89732e-09,-1.37569e-13,22706.3,-49.5769], Tmin=(1044.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C(O)C(=C)[CH]O[O](22891)',
    structure = SMILES('C=C(O)C(=C)[CH]O[O]'),
    E0 = (30.1067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.782041,0.0860902,-9.08788e-05,4.56122e-08,-8.57184e-12,3809.51,30.8839], Tmin=(100,'K'), Tmax=(1472.35,'K')), NASAPolynomial(coeffs=[23.8898,0.00844754,-9.6335e-07,2.50921e-12,3.98949e-15,-2304.99,-93.7807], Tmin=(1472.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.1067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C=C([CH]O[O])C(C)=O(22892)',
    structure = SMILES('C=C([CH]O[O])C(C)=O'),
    E0 = (17.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4951,0.0582329,-4.29841e-05,1.61575e-08,-2.56809e-12,2231.71,27.6204], Tmin=(100,'K'), Tmax=(1397.29,'K')), NASAPolynomial(coeffs=[10.4887,0.0324872,-1.53461e-05,2.97111e-09,-2.0883e-13,-281.654,-18.7774], Tmin=(1397.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C=C([O])C(=C)CO[O](2964)',
    structure = SMILES('C=C([O])C(=C)CO[O]'),
    E0 = (50.6164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,228.003,228.012,228.025],'cm^-1')),
        HinderedRotor(inertia=(0.306913,'amu*angstrom^2'), symmetry=1, barrier=(11.3259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306948,'amu*angstrom^2'), symmetry=1, barrier=(11.3257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30696,'amu*angstrom^2'), symmetry=1, barrier=(11.3259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283427,0.0825575,-0.000101641,6.61965e-08,-1.70641e-11,6220.89,27.9596], Tmin=(100,'K'), Tmax=(951.06,'K')), NASAPolynomial(coeffs=[14.1522,0.0242301,-9.65149e-06,1.71669e-09,-1.15273e-13,3582.78,-38.2537], Tmin=(951.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.6164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2][C]([O])C([CH2])[CH]O[O](22893)',
    structure = SMILES('[CH2][C]([O])C([CH2])[CH]O[O]'),
    E0 = (729.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.301586,0.105619,-0.000179581,1.56884e-07,-5.19527e-11,87904.2,35.3878], Tmin=(100,'K'), Tmax=(897.632,'K')), NASAPolynomial(coeffs=[10.447,0.0322412,-1.43831e-05,2.56926e-09,-1.67089e-13,87001.1,-9.58892], Tmin=(897.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH) + radical(CCsJOOH) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])[C]1CC1O[O](22894)',
    structure = SMILES('[CH2]C([O])[C]1CC1O[O]'),
    E0 = (496.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.947744,0.0568825,-2.88604e-05,-5.39189e-09,6.58298e-12,59882.1,33.6018], Tmin=(100,'K'), Tmax=(1018.28,'K')), NASAPolynomial(coeffs=[15.1843,0.0220081,-8.49585e-06,1.57648e-09,-1.12e-13,55891.5,-40.698], Tmin=(1018.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CJCOOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1OC[C]1[CH]O[O](22895)',
    structure = SMILES('[CH2]C1OC[C]1[CH]O[O]'),
    E0 = (488.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80052,0.0673091,-6.96307e-05,4.20755e-08,-1.01711e-11,58825.2,31.7848], Tmin=(100,'K'), Tmax=(1105.33,'K')), NASAPolynomial(coeffs=[10.7971,0.0265331,-8.05241e-06,1.17006e-09,-6.76311e-14,56896.3,-16.1728], Tmin=(1105.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOOH) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]O[CH][C]1CCC1[O](22896)',
    structure = SMILES('[O]O[CH][C]1CCC1[O]'),
    E0 = (482.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37879,0.0477046,-1.06669e-05,-1.61602e-08,8.26863e-12,58088.1,31.9331], Tmin=(100,'K'), Tmax=(1071.67,'K')), NASAPolynomial(coeffs=[12.1206,0.0280023,-1.16314e-05,2.195e-09,-1.55263e-13,54614.9,-26.0965], Tmin=(1071.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CJCOOH) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C([O])[C]1[CH]OOC1(22897)',
    structure = SMILES('[CH2]C([O])[C]1[CH]OOC1'),
    E0 = (465.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24798,0.0480979,5.19846e-06,-4.98175e-08,2.63155e-11,56134.7,30.9816], Tmin=(100,'K'), Tmax=(890.173,'K')), NASAPolynomial(coeffs=[15.432,0.0188747,-3.71505e-06,4.12555e-10,-2.33707e-14,52242,-43.4786], Tmin=(890.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CC(C)OJ) + radical(C2CJCOOH) + radical(CCsJOOC) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1OC1([CH2])[CH]O[O](22817)',
    structure = SMILES('[CH2]C1OC1([CH2])[CH]O[O]'),
    E0 = (492.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612488,0.103941,-0.000155308,1.1574e-07,-3.22919e-11,59427.6,29.2845], Tmin=(100,'K'), Tmax=(1032.93,'K')), NASAPolynomial(coeffs=[16.2046,0.0195233,-4.69979e-06,4.51635e-10,-1.20533e-14,56982.7,-47.4105], Tmin=(1032.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)OC) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1([CH]O[O])CC1[O](22752)',
    structure = SMILES('[CH2]C1([CH]O[O])CC1[O]'),
    E0 = (488.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.417308,0.0717682,-6.588e-05,3.04235e-08,-5.56516e-12,58895.5,29.5457], Tmin=(100,'K'), Tmax=(1320.05,'K')), NASAPolynomial(coeffs=[16.9733,0.0216002,-8.87304e-06,1.63317e-09,-1.12641e-13,54524.5,-54.9247], Tmin=(1320.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(ROOJ) + radical(CCsJOOH) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C([O])C1([CH2])[CH]OO1(22898)',
    structure = SMILES('[CH2]C([O])C1([CH2])[CH]OO1'),
    E0 = (553.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259569,0.0858316,-0.000107245,7.19152e-08,-1.93718e-11,66667.5,27.6907], Tmin=(100,'K'), Tmax=(904.301,'K')), NASAPolynomial(coeffs=[13.1004,0.0290335,-1.30332e-05,2.4612e-09,-1.7104e-13,64345.1,-32.9671], Tmin=(904.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(CC(C)OJ) + radical(CCsJOO) + radical(CJCOOH) + radical(CJCO)"""),
)

species(
    label = 'C=C([O])C(=C)[CH]O[O](19816)',
    structure = SMILES('C=C([O])C(=C)[CH]O[O]'),
    E0 = (167.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,365.271,365.782,366.4],'cm^-1')),
        HinderedRotor(inertia=(0.213009,'amu*angstrom^2'), symmetry=1, barrier=(20.2587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212055,'amu*angstrom^2'), symmetry=1, barrier=(20.2522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212784,'amu*angstrom^2'), symmetry=1, barrier=(20.2501,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4369.1,'J/mol'), sigma=(6.96718,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.44 K, Pc=29.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0353323,0.0777402,-8.35661e-05,4.40777e-08,-8.89542e-12,20349.8,29.6976], Tmin=(100,'K'), Tmax=(1286.73,'K')), NASAPolynomial(coeffs=[19.9352,0.012289,-3.33858e-06,4.75854e-10,-2.85611e-14,15489.4,-70.5996], Tmin=(1286.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C([O])C(=C)C=O(22899)',
    structure = SMILES('[CH2]C([O])C(=C)C=O'),
    E0 = (121.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,368.944,369.211],'cm^-1')),
        HinderedRotor(inertia=(0.00123805,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166126,'amu*angstrom^2'), symmetry=1, barrier=(15.999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166466,'amu*angstrom^2'), symmetry=1, barrier=(16.0057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2782,0.0598716,-5.04691e-05,2.10749e-08,-3.58269e-12,14748.7,26.18], Tmin=(100,'K'), Tmax=(1368.92,'K')), NASAPolynomial(coeffs=[13.262,0.0248546,-1.20991e-05,2.38857e-09,-1.70089e-13,11467.7,-35.3982], Tmin=(1368.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'C=C([CH]O[O])C=O(21078)',
    structure = SMILES('C=C([CH]O[O])C=O'),
    E0 = (62.9869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1712.55],'cm^-1')),
        HinderedRotor(inertia=(0.951473,'amu*angstrom^2'), symmetry=1, barrier=(21.8762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187831,'amu*angstrom^2'), symmetry=1, barrier=(21.8861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.954312,'amu*angstrom^2'), symmetry=1, barrier=(21.9415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80712,0.0497199,-4.02727e-05,1.59785e-08,-2.60191e-12,7652.85,23.5311], Tmin=(100,'K'), Tmax=(1405.63,'K')), NASAPolynomial(coeffs=[11.3369,0.0226005,-1.1332e-05,2.25213e-09,-1.60551e-13,4973.82,-25.6893], Tmin=(1405.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.9869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C=CC(=C)[CH]O[O](21014)',
    structure = SMILES('C=CC(=C)[CH]O[O]'),
    E0 = (244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,291.021,291.567],'cm^-1')),
        HinderedRotor(inertia=(0.393734,'amu*angstrom^2'), symmetry=1, barrier=(23.7609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395434,'amu*angstrom^2'), symmetry=1, barrier=(23.755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395156,'amu*angstrom^2'), symmetry=1, barrier=(23.7592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979578,0.0550171,-2.41496e-05,-1.47396e-08,1.14105e-11,29465.4,25.3096], Tmin=(100,'K'), Tmax=(965.846,'K')), NASAPolynomial(coeffs=[16.7953,0.0162363,-5.41757e-06,9.73081e-10,-7.0359e-14,25164,-56.8943], Tmin=(965.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C]=CO[O](16807)',
    structure = SMILES('C=[C][CH]O[O]'),
    E0 = (437.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,274.987],'cm^-1')),
        HinderedRotor(inertia=(0.170957,'amu*angstrom^2'), symmetry=1, barrier=(9.18865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171244,'amu*angstrom^2'), symmetry=1, barrier=(9.18809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2197,0.041819,-5.77934e-05,4.42117e-08,-1.35872e-11,52657.9,18.5046], Tmin=(100,'K'), Tmax=(795.764,'K')), NASAPolynomial(coeffs=[7.632,0.0146118,-6.50507e-06,1.24122e-09,-8.65948e-14,51796.6,-6.36981], Tmin=(795.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH][O](1556)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (367.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1667.93],'cm^-1')),
        HinderedRotor(inertia=(0.00517725,'amu*angstrom^2'), symmetry=1, barrier=(10.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22654,0.0212781,-3.59354e-05,3.91027e-08,-1.59281e-11,44278.5,12.2199], Tmin=(100,'K'), Tmax=(830.699,'K')), NASAPolynomial(coeffs=[2.17156,0.016018,-7.76586e-06,1.51127e-09,-1.05387e-13,44810.5,19.2612], Tmin=(830.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]([O])C(=C)[CH]O[O](22900)',
    structure = SMILES('[CH2]C([O])=C([CH2])[CH]O[O]'),
    E0 = (352.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,517.077,517.138],'cm^-1')),
        HinderedRotor(inertia=(0.172443,'amu*angstrom^2'), symmetry=1, barrier=(32.726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000630304,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052647,'amu*angstrom^2'), symmetry=1, barrier=(9.9926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.445701,'amu*angstrom^2'), symmetry=1, barrier=(84.5876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.2699,0.0749843,-7.39558e-05,3.1384e-08,-3.35158e-12,42521.5,31.5178], Tmin=(100,'K'), Tmax=(947.693,'K')), NASAPolynomial(coeffs=[18.3445,0.0149585,-4.6885e-06,7.64813e-10,-5.10911e-14,38365.4,-58.5637], Tmin=(947.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C([CH]O[O])C([CH2])[O](19807)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])[O]'),
    E0 = (648.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,895.683],'cm^-1')),
        HinderedRotor(inertia=(0.0956283,'amu*angstrom^2'), symmetry=1, barrier=(15.4523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672029,'amu*angstrom^2'), symmetry=1, barrier=(15.4513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210245,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271359,'amu*angstrom^2'), symmetry=1, barrier=(15.4522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202557,0.0774458,-8.31329e-05,4.38994e-08,-9.01416e-12,78118.2,34.4512], Tmin=(100,'K'), Tmax=(1194.96,'K')), NASAPolynomial(coeffs=[18.4095,0.0164997,-6.6286e-06,1.21751e-09,-8.4533e-14,73766.9,-56.6294], Tmin=(1194.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](O)C(=C)[CH]O[O](22901)',
    structure = SMILES('[CH2]C(O)=C([CH2])[CH]O[O]'),
    E0 = (214.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.738019,0.086426,-9.21004e-05,4.69378e-08,-8.9787e-12,25992.5,33.6399], Tmin=(100,'K'), Tmax=(1439.63,'K')), NASAPolynomial(coeffs=[23.4624,0.00919977,-1.23134e-06,3.98417e-11,2.09098e-15,20059.4,-88.3374], Tmin=(1439.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.569,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([CH]O[O])[C](C)[O](22902)',
    structure = SMILES('[CH2]C([CH]O[O])=C(C)[O]'),
    E0 = (193.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.17423,0.0759487,-7.57651e-05,3.82632e-08,-7.56384e-12,23412.2,31.1738], Tmin=(100,'K'), Tmax=(1238.51,'K')), NASAPolynomial(coeffs=[17.7848,0.0190715,-6.87841e-06,1.18242e-09,-7.87869e-14,19050.1,-57.5536], Tmin=(1238.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](CO[O])C(=C)[O](2965)',
    structure = SMILES('[CH2]C([O])=C([CH2])CO[O]'),
    E0 = (235.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,378.537,2470.92],'cm^-1')),
        HinderedRotor(inertia=(0.10079,'amu*angstrom^2'), symmetry=1, barrier=(10.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100788,'amu*angstrom^2'), symmetry=1, barrier=(10.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100791,'amu*angstrom^2'), symmetry=1, barrier=(10.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712185,'amu*angstrom^2'), symmetry=1, barrier=(72.4184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.245619,0.0838183,-0.000105923,7.12627e-08,-1.89768e-11,28407.6,31.012], Tmin=(100,'K'), Tmax=(921.727,'K')), NASAPolynomial(coeffs=[13.8602,0.0247353,-9.77183e-06,1.71844e-09,-1.14185e-13,25897.8,-33.5605], Tmin=(921.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(CO[O])C([CH2])[O](19309)',
    structure = SMILES('[CH]=C(CO[O])C([CH2])[O]'),
    E0 = (531.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,348.945,348.946],'cm^-1')),
        HinderedRotor(inertia=(0.00138453,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116446,'amu*angstrom^2'), symmetry=1, barrier=(10.0617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11645,'amu*angstrom^2'), symmetry=1, barrier=(10.0617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116446,'amu*angstrom^2'), symmetry=1, barrier=(10.0617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4372.93,'J/mol'), sigma=(7.19133,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=683.04 K, Pc=26.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355041,0.0841215,-0.000107279,7.34023e-08,-2.0171e-11,63996.7,33.3171], Tmin=(100,'K'), Tmax=(886.899,'K')), NASAPolynomial(coeffs=[12.7302,0.0283088,-1.2885e-05,2.44852e-09,-1.70658e-13,61801.6,-24.9008], Tmin=(886.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C([CH2])O(20473)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])O'),
    E0 = (417.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,370.456],'cm^-1')),
        HinderedRotor(inertia=(0.0110901,'amu*angstrom^2'), symmetry=1, barrier=(20.0921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873887,'amu*angstrom^2'), symmetry=1, barrier=(20.0924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00408916,'amu*angstrom^2'), symmetry=1, barrier=(20.0923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873899,'amu*angstrom^2'), symmetry=1, barrier=(20.0927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873886,'amu*angstrom^2'), symmetry=1, barrier=(20.0924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0749316,0.0810327,-8.76603e-05,4.72219e-08,-9.92549e-12,50416.1,35.1008], Tmin=(100,'K'), Tmax=(1166.41,'K')), NASAPolynomial(coeffs=[18.2735,0.0186229,-7.40009e-06,1.34802e-09,-9.30364e-14,46170.8,-55.4975], Tmin=(1166.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C(C)[O](20474)',
    structure = SMILES('[CH]=C([CH]O[O])C(C)[O]'),
    E0 = (436.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,367.44,890.006],'cm^-1')),
        HinderedRotor(inertia=(0.137313,'amu*angstrom^2'), symmetry=1, barrier=(3.15711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896564,'amu*angstrom^2'), symmetry=1, barrier=(20.6138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.896538,'amu*angstrom^2'), symmetry=1, barrier=(20.6132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0366694,'amu*angstrom^2'), symmetry=1, barrier=(20.6155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.169993,0.0743386,-7.01906e-05,3.28696e-08,-6.01977e-12,52674.2,33.9403], Tmin=(100,'K'), Tmax=(1331.59,'K')), NASAPolynomial(coeffs=[18.824,0.0183026,-7.06701e-06,1.26614e-09,-8.62927e-14,47706.3,-61.3964], Tmin=(1331.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]([O])C(=C)[CH]OO(22903)',
    structure = SMILES('[CH2]C([O])=C([CH2])[CH]OO'),
    E0 = (200.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.434815,0.0836644,-8.65047e-05,4.36325e-08,-8.42025e-12,24270.4,33.2944], Tmin=(100,'K'), Tmax=(1352.41,'K')), NASAPolynomial(coeffs=[21.9432,0.0132827,-3.78963e-06,5.64993e-10,-3.50462e-14,18601.1,-80.004], Tmin=(1352.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C([CH]OO)C([CH2])[O](20471)',
    structure = SMILES('[CH]=C([CH]OO)C([CH2])[O]'),
    E0 = (496.327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,339.742],'cm^-1')),
        HinderedRotor(inertia=(1.11562,'amu*angstrom^2'), symmetry=1, barrier=(26.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02149,'amu*angstrom^2'), symmetry=1, barrier=(23.486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0294512,'amu*angstrom^2'), symmetry=1, barrier=(23.5594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0294063,'amu*angstrom^2'), symmetry=1, barrier=(23.5274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0293689,'amu*angstrom^2'), symmetry=1, barrier=(23.5072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.218476,0.082812,-8.42355e-05,4.15046e-08,-7.92122e-12,59854.6,35.208], Tmin=(100,'K'), Tmax=(1285.51,'K')), NASAPolynomial(coeffs=[21.1845,0.0162146,-6.5266e-06,1.20493e-09,-8.39511e-14,54351.8,-73.4247], Tmin=(1285.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CCJO) + radical(CJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C[C]=CO[O](21021)',
    structure = SMILES('[CH2]C([O])C[C]=CO[O]'),
    E0 = (534.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,339.822,339.822,339.822],'cm^-1')),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104059,'amu*angstrom^2'), symmetry=1, barrier=(8.52727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4429.24,'J/mol'), sigma=(7.2345,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.84 K, Pc=26.54 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28743,0.0846749,-0.00010982,7.56961e-08,-2.07614e-11,64443.7,33.5907], Tmin=(100,'K'), Tmax=(892.541,'K')), NASAPolynomial(coeffs=[13.3715,0.0260365,-1.12712e-05,2.08583e-09,-1.42994e-13,62108.1,-28.045], Tmin=(892.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])[C]=CO[O](22904)',
    structure = SMILES('[CH2]C([O])[C]=CO[O]'),
    E0 = (558.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,314.9,314.93],'cm^-1')),
        HinderedRotor(inertia=(0.133849,'amu*angstrom^2'), symmetry=1, barrier=(9.42344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133969,'amu*angstrom^2'), symmetry=1, barrier=(9.42188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13406,'amu*angstrom^2'), symmetry=1, barrier=(9.42349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.973005,0.0668437,-8.43773e-05,5.38556e-08,-1.34401e-11,67241.2,29.79], Tmin=(100,'K'), Tmax=(985.519,'K')), NASAPolynomial(coeffs=[13.5725,0.015704,-6.539e-06,1.19991e-09,-8.24242e-14,64757.8,-30.8114], Tmin=(985.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1OCC1=CO[O](22905)',
    structure = SMILES('[CH2]C1OCC1=CO[O]'),
    E0 = (245.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732015,0.0613261,-3.63762e-05,-1.22659e-09,5.97303e-12,29649.1,29.5427], Tmin=(100,'K'), Tmax=(1004.36,'K')), NASAPolynomial(coeffs=[16.6229,0.0200116,-7.49039e-06,1.38291e-09,-9.8682e-14,25348.8,-52.7082], Tmin=(1004.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]OC=C1CCC1[O](22704)',
    structure = SMILES('[O]OC=C1CCC1[O]'),
    E0 = (241.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02749,0.053912,-1.89661e-05,-1.55143e-08,9.97337e-12,29135.7,28.6713], Tmin=(100,'K'), Tmax=(1016.26,'K')), NASAPolynomial(coeffs=[15.2999,0.0223728,-8.77842e-06,1.65771e-09,-1.19323e-13,24962.6,-46.6743], Tmin=(1016.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([O])C1=COOC1(19173)',
    structure = SMILES('[CH2]C([O])C1=COOC1'),
    E0 = (213.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771971,0.0564396,-1.41182e-05,-2.86671e-08,1.64965e-11,25758.5,30.4679], Tmin=(100,'K'), Tmax=(980.519,'K')), NASAPolynomial(coeffs=[18.6067,0.0175085,-6.30705e-06,1.20443e-09,-9.01362e-14,20635,-63.5149], Tmin=(980.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1=COOOC1[CH2](22906)',
    structure = SMILES('[CH2]C1OOO[CH]C1=C'),
    E0 = (281.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39242,0.0230889,0.000110415,-1.68788e-07,6.71135e-11,34037.4,26.9274], Tmin=(100,'K'), Tmax=(971.341,'K')), NASAPolynomial(coeffs=[25.2363,0.0110456,-4.01869e-06,1.05713e-09,-1.00236e-13,25341.3,-108.332], Tmin=(971.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(C=CCJO) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1=COOCC1[O](22907)',
    structure = SMILES('C=C1[CH]OOCC1[O]'),
    E0 = (102.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27566,0.0323087,7.40122e-05,-1.22325e-07,4.8343e-11,12457.6,23.7501], Tmin=(100,'K'), Tmax=(996.033,'K')), NASAPolynomial(coeffs=[21.6978,0.0187243,-8.58315e-06,1.93388e-09,-1.58048e-13,4994.98,-91.7334], Tmin=(996.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(CC(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C([O])C1=COC1(20462)',
    structure = SMILES('[CH2]C([O])C1=COC1'),
    E0 = (212.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930969,0.0446241,2.98795e-05,-8.83898e-08,4.19383e-11,25660.9,24.4389], Tmin=(100,'K'), Tmax=(929.754,'K')), NASAPolynomial(coeffs=[24.1811,0.00329574,1.85542e-06,-3.91923e-10,1.82595e-14,18800.4,-99.68], Tmin=(929.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1=COOC1[CH2](22908)',
    structure = SMILES('[CH2]C1OO[CH]C1=C'),
    E0 = (267.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96782,0.0171298,9.89974e-05,-1.47863e-07,5.92474e-11,32261.3,25.1591], Tmin=(100,'K'), Tmax=(951.206,'K')), NASAPolynomial(coeffs=[19.8555,0.0117532,-2.66524e-06,5.8321e-10,-5.62777e-14,25698.6,-76.8531], Tmin=(951.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJO) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1=COCC1[O](22909)',
    structure = SMILES('[CH2]C1=COCC1[O]'),
    E0 = (41.4745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48715,0.0298838,6.94775e-05,-1.30449e-07,5.83059e-11,5102.5,18.78], Tmin=(100,'K'), Tmax=(905.735,'K')), NASAPolynomial(coeffs=[22.9301,0.00213048,4.57101e-06,-1.06893e-09,7.01151e-14,-1527.79,-97.7058], Tmin=(905.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.4745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=O)C(C)=CO[O](22910)',
    structure = SMILES('C=C([O])C(C)=CO[O]'),
    E0 = (47.9494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00436091,0.0843218,-0.000101429,6.20122e-08,-1.4706e-11,5913.95,28.3405], Tmin=(100,'K'), Tmax=(1042.13,'K')), NASAPolynomial(coeffs=[17.2771,0.0180239,-6.0023e-06,9.65872e-10,-6.13094e-14,2313.89,-55.7031], Tmin=(1042.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.9494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(=O)C([CH2])=COO(22911)',
    structure = SMILES('C=C([O])C(=C)[CH]OO'),
    E0 = (15.9069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0175901,0.0780406,-6.75007e-05,2.02664e-08,9.47295e-13,2067.06,28.8736], Tmin=(100,'K'), Tmax=(974.175,'K')), NASAPolynomial(coeffs=[20.6089,0.0154563,-5.17823e-06,9.13465e-10,-6.47941e-14,-3000.81,-75.4824], Tmin=(974.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.9069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C]1C([CH2])OC1O[O](22912)',
    structure = SMILES('[CH2][C]1C([CH2])OC1O[O]'),
    E0 = (464.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528925,0.0654365,-6.11354e-05,3.1796e-08,-6.40209e-12,55966.8,34.1524], Tmin=(100,'K'), Tmax=(1401.42,'K')), NASAPolynomial(coeffs=[12.8591,0.0214747,-4.69568e-06,4.82487e-10,-1.95776e-14,53371.9,-26.4225], Tmin=(1401.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(C2CJCOOH) + radical(CJC(C)OC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1C([O])CC1O[O](22913)',
    structure = SMILES('[CH2][C]1C([O])CC1O[O]'),
    E0 = (485.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4338,0.0440355,5.75417e-06,-3.92134e-08,1.82817e-11,58552.7,32.3255], Tmin=(100,'K'), Tmax=(975.359,'K')), NASAPolynomial(coeffs=[13.5429,0.0238566,-8.55155e-06,1.55408e-09,-1.10759e-13,54788.2,-32.9803], Tmin=(975.359,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C([CH2])[O])C1OO1(22914)',
    structure = SMILES('[CH2][C](C([CH2])[O])C1OO1'),
    E0 = (525.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.857797,0.0568145,-1.69475e-05,-2.94822e-08,1.91789e-11,63369.3,35.6862], Tmin=(100,'K'), Tmax=(909.241,'K')), NASAPolynomial(coeffs=[17.7682,0.0154555,-3.21442e-06,4.07136e-10,-2.61898e-14,58928.7,-51.797], Tmin=(909.241,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(dioxirane) + radical(CC(C)OJ) + radical(C2CJCOOH) + radical(CJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=[C]OO)C([CH2])[O](22915)',
    structure = SMILES('[CH2]C(=[C]OO)C([CH2])[O]'),
    E0 = (520.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00705025,'amu*angstrom^2'), symmetry=1, barrier=(4.54811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.677851,'amu*angstrom^2'), symmetry=1, barrier=(15.5851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0241468,'amu*angstrom^2'), symmetry=1, barrier=(15.5834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0241273,'amu*angstrom^2'), symmetry=1, barrier=(15.586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6491,'amu*angstrom^2'), symmetry=1, barrier=(37.9161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277306,0.0830355,-9.61826e-05,5.73638e-08,-1.36121e-11,62736.3,35.2247], Tmin=(100,'K'), Tmax=(1025.38,'K')), NASAPolynomial(coeffs=[15.1566,0.0249911,-1.12708e-05,2.15688e-09,-1.51907e-13,59684.9,-36.9322], Tmin=(1025.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]([O])C(C)=CO[O](22916)',
    structure = SMILES('[CH2]C([O])=C(C)[CH]O[O]'),
    E0 = (200.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.237618,0.0778677,-8.38906e-05,4.65536e-08,-1.01361e-11,24299.4,31.3527], Tmin=(100,'K'), Tmax=(1126.51,'K')), NASAPolynomial(coeffs=[16.4286,0.0203777,-7.34108e-06,1.25234e-09,-8.28406e-14,20651.4,-48.6886], Tmin=(1126.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])C(C)=[C]O[O](22917)',
    structure = SMILES('[CH2]C([O])C(C)=[C]O[O]'),
    E0 = (521.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.432721,0.083253,-0.000114274,8.68372e-08,-2.66165e-11,62787.9,35.1438], Tmin=(100,'K'), Tmax=(797.097,'K')), NASAPolynomial(coeffs=[11.1038,0.0296991,-1.34877e-05,2.53639e-09,-1.74637e-13,61086.9,-13.9174], Tmin=(797.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C([CH2])O(22918)',
    structure = SMILES('[CH2]C(=[C]O[O])C([CH2])O'),
    E0 = (442.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,401.028],'cm^-1')),
        HinderedRotor(inertia=(0.0792732,'amu*angstrom^2'), symmetry=1, barrier=(9.00209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0796141,'amu*angstrom^2'), symmetry=1, barrier=(8.99176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0792691,'amu*angstrom^2'), symmetry=1, barrier=(8.99085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0794126,'amu*angstrom^2'), symmetry=1, barrier=(9.00309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00104962,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.303274,0.0843171,-0.000109864,7.57452e-08,-2.07523e-11,53309.5,36.0829], Tmin=(100,'K'), Tmax=(893.839,'K')), NASAPolynomial(coeffs=[13.4611,0.0254359,-1.10544e-05,2.04985e-09,-1.40701e-13,50957.3,-25.9197], Tmin=(893.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CJCO) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=[C]O[O])C(C)[O](22919)',
    structure = SMILES('[CH2]C(=[C]O[O])C(C)[O]'),
    E0 = (460.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,396.915,397.02],'cm^-1')),
        HinderedRotor(inertia=(0.00106978,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0836922,'amu*angstrom^2'), symmetry=1, barrier=(9.35948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0837313,'amu*angstrom^2'), symmetry=1, barrier=(9.35808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.083713,'amu*angstrom^2'), symmetry=1, barrier=(9.35896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690232,0.0743004,-8.1341e-05,4.78525e-08,-1.1399e-11,55554.8,33.8675], Tmin=(100,'K'), Tmax=(1013.93,'K')), NASAPolynomial(coeffs=[12.5713,0.0274292,-1.20003e-05,2.26061e-09,-1.57645e-13,53145.4,-23.6161], Tmin=(1013.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(ROOJ) + radical(Allyl_P) + radical(C=CJO)"""),
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
    E0 = (401.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (558.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (612.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (817.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (631.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (941.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (676.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1016.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (884.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (849.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (409.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (406.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (409.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (409.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (424.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (424.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (464.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (752.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (632.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (527.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (531.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (465.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (496.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (488.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (555.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (401.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (401.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (464.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (502.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (502.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (642.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (805.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (564.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (860.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (515.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (542.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (577.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (676.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (462.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (481.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (494.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (529.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (704.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (995.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (409.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (409.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (408.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (409.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (410.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (484.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (510.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (463.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (464.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (426.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (527.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (531.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (525.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (662.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (603.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (656.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (486.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (505.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['vinoxy(1351)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C([CH]O[O])C[CH][O](20382)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([O])C([C]=C)O[O](21020)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH2]C([O])C(=C)[CH][O](22887)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O2(2)', '[CH]=C([CH2])C([CH2])[O](15788)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2][CH]C(=C)[CH]O[O](21398)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', 'C=C([CH][O])[CH]O[O](22888)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]O[O](21387)', '[CH2]C([O])[C]=C(2851)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2]C(=[C]O[O])C([CH2])[O](22522)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([O])C(=C)[CH]O[O](21139)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C([O])C(=C)C1OO1(22889)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C([CH]O[O])C1CO1(22890)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1OC(O[O])C1=C(22791)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C1C([O])CC1O[O](22719)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C(O)C(=C)[CH]O[O](22891)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C([CH]O[O])C(C)=O(22892)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C([O])C(=C)CO[O](2964)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]([O])C([CH2])[CH]O[O](22893)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C([O])[C]1CC1O[O](22894)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1OC[C]1[CH]O[O](22895)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[O]O[CH][C]1CCC1[O](22896)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C([O])[C]1[CH]OOC1(22897)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.80789e+10,'s^-1'), n=0.280222, Ea=(64.5731,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 61.3 to 64.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1OC1([CH2])[CH]O[O](22817)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(95.4145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1([CH]O[O])CC1[O](22752)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(87.3284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 86.9 to 87.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C([O])C1([CH2])[CH]OO1(22898)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(153.943,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;doublebond_intra_2H;radadd_intra_O] + [R5;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'C=C([O])C(=C)[CH]O[O](19816)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(4)', '[CH2]C([O])C(=C)C=O(22899)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(36.4151,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 32.2 to 36.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(T)(20)', 'C=C([CH]O[O])C=O(21078)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.5e+06,'cm^3/(mol*s)'), n=2.16, Ea=(19.2464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O(4)', 'C=CC(=C)[CH]O[O](21014)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.0532,'m^3/(mol*s)'), n=2.065, Ea=(15.9938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;YJ] for rate rule [Cds-CdH_Cds-HH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['vinoxy(1351)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][O](1556)', 'C=C=CO[O](16806)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00267131,'m^3/(mol*s)'), n=2.569, Ea=(24.4155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][O](1556)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2][C]([O])C(=C)[CH]O[O](22900)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH]=C([CH]O[O])C([CH2])[O](19807)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C](O)C(=C)[CH]O[O](22901)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['C=C([CH]O[O])[C](C)[O](22902)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C](CO[O])C(=C)[O](2965)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.00351592,'s^-1'), n=4.63833, Ea=(175.937,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeO;XH_out] + [R3H_SS_2Cd;C_rad_out_1H;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C(CO[O])C([CH2])[O](19309)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C([CH]O[O])C([CH2])O(20473)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C([CH]O[O])C(C)[O](20474)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C]([O])C(=C)[CH]OO(22903)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(36450.5,'s^-1'), n=1.97932, Ea=(93.5766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;XH_out] for rate rule [R5HJ_2;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C([CH]OO)C([CH2])[O](20471)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])C[C]=CO[O](21021)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CH2(T)(20)', '[CH2]C([O])[C]=CO[O](22904)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1OCC1=CO[O](22905)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[O]OC=C1CCC1[O](22704)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C([O])C1=COOC1(19173)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1=COOOC1[CH2](22906)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSDSS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C1=COOCC1[O](22907)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSDSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['O(4)', '[CH2]C([O])C1=COC1(20462)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['O(4)', '[CH2]C1=COOC1[CH2](22908)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(109.188,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO] for rate rule [R4OO_SSD;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation
Ea raised from 107.8 to 109.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['O(4)', '[CH2]C1=COCC1[O](22909)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(5.13e+10,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;C_pri_rad_intra;OO] for rate rule [R4OO_SSD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C(=O)C(C)=CO[O](22910)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C(=O)C([CH2])=COO(22911)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C]1C([CH2])OC1O[O](22912)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C]1C([O])CC1O[O](22913)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(130.445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C](C([CH2])[O])C1OO1(22914)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(124.611,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 122.4 to 124.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2]C(=[C]OO)C([CH2])[O](22915)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2][C]([O])C(C)=CO[O](22916)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(5.4947e+07,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    products = ['[CH2]C([O])C(C)=[C]O[O](22917)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C(=[C]O[O])C([CH2])O(22918)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C(=[C]O[O])C(C)[O](22919)'],
    products = ['[CH2]C([O])C(=C)[CH]O[O](20472)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4694',
    isomers = [
        '[CH2]C([O])C(=C)[CH]O[O](20472)',
    ],
    reactants = [
        ('vinoxy(1351)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4694',
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


species(
    label = '[CH]=C(CO[O])C([CH2])[CH2](19280)',
    structure = SMILES('[CH]=C(CO[O])C([CH2])[CH2]'),
    E0 = (644.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,1807.41],'cm^-1')),
        HinderedRotor(inertia=(0.207387,'amu*angstrom^2'), symmetry=1, barrier=(4.76824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207368,'amu*angstrom^2'), symmetry=1, barrier=(4.7678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207726,'amu*angstrom^2'), symmetry=1, barrier=(4.77603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207265,'amu*angstrom^2'), symmetry=1, barrier=(4.76542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84817,'amu*angstrom^2'), symmetry=1, barrier=(65.4851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153661,0.0872739,-0.000110947,8.07318e-08,-2.35617e-11,77650.1,35.1899], Tmin=(100,'K'), Tmax=(901.134,'K')), NASAPolynomial(coeffs=[11.3133,0.0329827,-1.26601e-05,2.16282e-09,-1.39884e-13,75831.9,-16.416], Tmin=(901.134,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C#CCO[O](16808)',
    structure = SMILES('C#CCO[O]'),
    E0 = (242.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.595127,'amu*angstrom^2'), symmetry=1, barrier=(13.6831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74124,'amu*angstrom^2'), symmetry=1, barrier=(63.0265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3504.08,'J/mol'), sigma=(5.7666,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=547.33 K, Pc=41.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29019,0.0394598,-5.33344e-05,4.08815e-08,-1.22872e-11,29219.5,16.1946], Tmin=(100,'K'), Tmax=(942.252,'K')), NASAPolynomial(coeffs=[6.8169,0.0148149,-5.46001e-06,8.95137e-10,-5.57485e-14,28607.4,-4.09622], Tmin=(942.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = 'allyl(82)',
    structure = SMILES('[CH2]C=C'),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29615,0.00579223,4.3392e-05,-5.9989e-08,2.33815e-11,18908.2,9.01994], Tmin=(100,'K'), Tmax=(942.181,'K')), NASAPolynomial(coeffs=[8.06863,0.0101837,-2.84795e-06,5.0088e-10,-3.79628e-14,16914.7,-19.5272], Tmin=(942.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C(C[CH][CH2])CO[O](19278)',
    structure = SMILES('[CH]=C(C[CH][CH2])CO[O]'),
    E0 = (637.806,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,1878.74],'cm^-1')),
        HinderedRotor(inertia=(0.144323,'amu*angstrom^2'), symmetry=1, barrier=(3.31826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1443,'amu*angstrom^2'), symmetry=1, barrier=(3.31774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144299,'amu*angstrom^2'), symmetry=1, barrier=(3.31772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144258,'amu*angstrom^2'), symmetry=1, barrier=(3.31677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144384,'amu*angstrom^2'), symmetry=1, barrier=(3.31968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428907,0.0866649,-0.000124421,1.08273e-07,-3.77952e-11,76831.1,36.0976], Tmin=(100,'K'), Tmax=(838.746,'K')), NASAPolynomial(coeffs=[6.40031,0.0431533,-1.97194e-05,3.68207e-09,-2.50833e-13,76358.2,11.4915], Tmin=(838.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(RCCJC) + radical(RCCJ) + radical(Cds_P)"""),
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
    label = '[CH]=C(C[O])C([CH2])[CH2](20009)',
    structure = SMILES('[CH]=C(C[O])C([CH2])[CH2]'),
    E0 = (646.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,180,1650.99],'cm^-1')),
        HinderedRotor(inertia=(0.0328602,'amu*angstrom^2'), symmetry=1, barrier=(63.6164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00302086,'amu*angstrom^2'), symmetry=1, barrier=(5.83579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253886,'amu*angstrom^2'), symmetry=1, barrier=(5.83733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253864,'amu*angstrom^2'), symmetry=1, barrier=(5.83684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959135,0.0708025,-8.59387e-05,6.61233e-08,-2.09902e-11,77884,31.6397], Tmin=(100,'K'), Tmax=(871.463,'K')), NASAPolynomial(coeffs=[6.97588,0.0363121,-1.47412e-05,2.60669e-09,-1.72489e-13,77096.3,4.93767], Tmin=(871.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(CCOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]=C([CH][CH2])CO[O](19537)',
    structure = SMILES('[CH]C(=C[CH2])CO[O]'),
    E0 = (523.601,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,445.995,446.365,446.466,446.529],'cm^-1')),
        HinderedRotor(inertia=(0.379141,'amu*angstrom^2'), symmetry=1, barrier=(53.5896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379362,'amu*angstrom^2'), symmetry=1, barrier=(53.5883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379433,'amu*angstrom^2'), symmetry=1, barrier=(53.5917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379449,'amu*angstrom^2'), symmetry=1, barrier=(53.5871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17643,0.065264,-6.05829e-05,3.36983e-08,-8.16429e-12,63073.6,26.5396], Tmin=(100,'K'), Tmax=(963.074,'K')), NASAPolynomial(coeffs=[8.01698,0.0368527,-1.63322e-05,3.0668e-09,-2.12821e-13,61756,-6.20478], Tmin=(963.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH2])C(=[CH])CO[O](20010)',
    structure = SMILES('[CH]C([CH2])C(=[CH])CO[O]'),
    E0 = (887.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320445,0.0855125,-0.00010521,6.82187e-08,-1.61579e-11,106885,34.5276], Tmin=(100,'K'), Tmax=(714.62,'K')), NASAPolynomial(coeffs=[11.8437,0.0312983,-1.30035e-05,2.34025e-09,-1.57303e-13,104975,-19.031], Tmin=(714.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(887.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=C(CO[O])C([CH2])[CH2](20011)',
    structure = SMILES('[C]=C(CO[O])C([CH2])[CH2]'),
    E0 = (955.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,845.21],'cm^-1')),
        HinderedRotor(inertia=(0.154427,'amu*angstrom^2'), symmetry=1, barrier=(3.55058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154322,'amu*angstrom^2'), symmetry=1, barrier=(3.54818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15428,'amu*angstrom^2'), symmetry=1, barrier=(3.5472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154376,'amu*angstrom^2'), symmetry=1, barrier=(3.54941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154192,'amu*angstrom^2'), symmetry=1, barrier=(3.54518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.187133,0.0901772,-0.000131691,1.07123e-07,-3.4019e-11,115051,35.0442], Tmin=(100,'K'), Tmax=(903.178,'K')), NASAPolynomial(coeffs=[9.76582,0.0334095,-1.35858e-05,2.35896e-09,-1.52524e-13,113906,-6.9528], Tmin=(903.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(955.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C(CO[O])C1CC1(20012)',
    structure = SMILES('[CH]=C(CO[O])C1CC1'),
    E0 = (392.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.642057,0.0639111,-3.57661e-05,9.82389e-10,4.13462e-12,47364.9,30.1201], Tmin=(100,'K'), Tmax=(1045.8,'K')), NASAPolynomial(coeffs=[14.8319,0.0280929,-1.08622e-05,1.98114e-09,-1.37936e-13,43387.7,-43.798], Tmin=(1045.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1CC=C1CO[O](20013)',
    structure = SMILES('[CH2]C1CC=C1CO[O]'),
    E0 = (349.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865264,0.0594352,-2.55829e-05,-8.43902e-09,7.5333e-12,42110.3,29.6715], Tmin=(100,'K'), Tmax=(997.999,'K')), NASAPolynomial(coeffs=[13.3934,0.0293935,-1.0748e-05,1.90371e-09,-1.3083e-13,38605.1,-35.7774], Tmin=(997.999,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C1=COOC1(19163)',
    structure = SMILES('[CH2]C([CH2])C1=COOC1'),
    E0 = (322.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.842246,0.0562368,-5.36599e-06,-3.83661e-08,2.07546e-11,38897.1,31.3766], Tmin=(100,'K'), Tmax=(933.691,'K')), NASAPolynomial(coeffs=[16.1665,0.0239578,-7.12113e-06,1.16678e-09,-8.01045e-14,34580.9,-49.2925], Tmin=(933.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1COOCC1[CH2](18244)',
    structure = SMILES('[CH]=C1COOCC1[CH2]'),
    E0 = (348.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903303,0.0500668,1.93584e-05,-5.81002e-08,2.4156e-11,41998.3,25.5994], Tmin=(100,'K'), Tmax=(1030.77,'K')), NASAPolynomial(coeffs=[16.8471,0.0303324,-1.32429e-05,2.64466e-09,-1.95864e-13,36472.9,-62.6615], Tmin=(1030.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH2])C1=COC1(20014)',
    structure = SMILES('[CH2]C([CH2])C1=COC1'),
    E0 = (321.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.992297,0.044484,3.86562e-05,-9.85653e-08,4.66423e-11,38799.9,25.3822], Tmin=(100,'K'), Tmax=(904.313,'K')), NASAPolynomial(coeffs=[21.9686,0.0093658,1.2568e-06,-4.79905e-10,3.24313e-14,32648.2,-86.7437], Tmin=(904.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1COCC1[CH2](18246)',
    structure = SMILES('[CH]=C1COCC1[CH2]'),
    E0 = (318.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82343,0.0316774,4.59538e-05,-8.22063e-08,3.41256e-11,38445.9,25.0247], Tmin=(100,'K'), Tmax=(944.14,'K')), NASAPolynomial(coeffs=[13.1069,0.0246237,-7.58162e-06,1.31062e-09,-9.38262e-14,34499,-38.3814], Tmin=(944.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(CO[O])C(=C)C(19572)',
    structure = SMILES('[CH]=C(CO[O])C(=C)C'),
    E0 = (336.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.570744,'amu*angstrom^2'), symmetry=1, barrier=(13.1225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570766,'amu*angstrom^2'), symmetry=1, barrier=(13.123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.57065,'amu*angstrom^2'), symmetry=1, barrier=(13.1204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570533,'amu*angstrom^2'), symmetry=1, barrier=(13.1177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229796,0.0850468,-9.48887e-05,5.8264e-08,-1.45432e-11,40570.7,28.4167], Tmin=(100,'K'), Tmax=(968.949,'K')), NASAPolynomial(coeffs=[12.8916,0.0327764,-1.39706e-05,2.58983e-09,-1.78655e-13,38117,-32.27], Tmin=(968.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=C)C(=C)CO[O](16266)',
    structure = SMILES('[CH2]C(=C)C(=C)CO[O]'),
    E0 = (240.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3000,3100,440,815,1455,1000,321.124,325.681],'cm^-1')),
        HinderedRotor(inertia=(0.17998,'amu*angstrom^2'), symmetry=1, barrier=(13.4205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180304,'amu*angstrom^2'), symmetry=1, barrier=(13.4089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507294,'amu*angstrom^2'), symmetry=1, barrier=(37.9631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180798,'amu*angstrom^2'), symmetry=1, barrier=(13.4096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0817636,0.0817437,-7.95512e-05,4.06894e-08,-8.32661e-12,29083.9,28.1239], Tmin=(100,'K'), Tmax=(1181.93,'K')), NASAPolynomial(coeffs=[16.0718,0.0276283,-1.08726e-05,1.95115e-09,-1.32717e-13,25304.1,-51.6917], Tmin=(1181.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(COO)C([CH2])=C(20015)',
    structure = SMILES('[CH]=C(COO)C([CH2])=C'),
    E0 = (335.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.236216,0.0886095,-9.05465e-05,4.7377e-08,-9.82054e-12,40532.3,28.9104], Tmin=(100,'K'), Tmax=(1172.45,'K')), NASAPolynomial(coeffs=[18.0642,0.0261727,-1.06642e-05,1.95364e-09,-1.34678e-13,36241.2,-62.2896], Tmin=(1172.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(CO[O])C([CH2])=C(19489)',
    structure = SMILES('[CH]=C(CO[O])C([CH2])=C'),
    E0 = (487.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,222.496],'cm^-1')),
        HinderedRotor(inertia=(0.424134,'amu*angstrom^2'), symmetry=1, barrier=(14.8406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423787,'amu*angstrom^2'), symmetry=1, barrier=(14.8301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419331,'amu*angstrom^2'), symmetry=1, barrier=(14.8292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.37592,'amu*angstrom^2'), symmetry=1, barrier=(116.708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102827,0.0841266,-9.21653e-05,5.2826e-08,-1.20325e-11,58799.6,28.4536], Tmin=(100,'K'), Tmax=(1069.76,'K')), NASAPolynomial(coeffs=[15.6908,0.0258405,-1.04367e-05,1.89293e-09,-1.29462e-13,55464.6,-47.8002], Tmin=(1069.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)CO[O](19302)',
    structure = SMILES('[CH]=C(C=C)CO[O]'),
    E0 = (373.801,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.668807,'amu*angstrom^2'), symmetry=1, barrier=(15.3772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668639,'amu*angstrom^2'), symmetry=1, barrier=(15.3733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.669376,'amu*angstrom^2'), symmetry=1, barrier=(15.3903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.929731,0.0668241,-7.09229e-05,4.01811e-08,-9.13985e-12,45069,25.2974], Tmin=(100,'K'), Tmax=(1065.74,'K')), NASAPolynomial(coeffs=[12.5691,0.023139,-9.4382e-06,1.72027e-09,-1.17892e-13,42588.1,-31.5972], Tmin=(1065.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CO[O](16805)',
    structure = SMILES('[CH]=[C]CO[O]'),
    E0 = (567.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.348758,'amu*angstrom^2'), symmetry=1, barrier=(8.01863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349832,'amu*angstrom^2'), symmetry=1, barrier=(8.04334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92348,0.0566508,-0.000115722,1.14457e-07,-4.10713e-11,68272.1,19.3685], Tmin=(100,'K'), Tmax=(892.364,'K')), NASAPolynomial(coeffs=[4.05419,0.0203884,-9.86661e-06,1.83106e-09,-1.20972e-13,68955.4,15.2908], Tmin=(892.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH][CH2](497)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.000863049,'amu*angstrom^2'), symmetry=1, barrier=(2.40754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000864751,'amu*angstrom^2'), symmetry=1, barrier=(2.41365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42318,0.0135015,1.25385e-06,-4.52647e-09,1.21417e-12,58330.2,15.4092], Tmin=(100,'K'), Tmax=(1578.37,'K')), NASAPolynomial(coeffs=[5.16145,0.0152713,-6.29635e-06,1.14122e-09,-7.61383e-14,57012.3,3.79304], Tmin=(1578.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]O[O](46)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (205.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2028.66,2028.67],'cm^-1')),
        HinderedRotor(inertia=(0.221535,'amu*angstrom^2'), symmetry=1, barrier=(13.5753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43223,0.0119348,-9.20111e-06,4.32805e-09,-8.54519e-13,24777.3,10.782], Tmin=(100,'K'), Tmax=(1197.14,'K')), NASAPolynomial(coeffs=[5.14053,0.00622677,-2.04895e-06,3.45073e-10,-2.27385e-14,24368.2,2.23311], Tmin=(1197.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = 'C#CC([CH2])[CH2](9851)',
    structure = SMILES('C#CC([CH2])[CH2]'),
    E0 = (526.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.218675,'amu*angstrom^2'), symmetry=1, barrier=(75.8207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.2973,'amu*angstrom^2'), symmetry=1, barrier=(75.8114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0119429,'amu*angstrom^2'), symmetry=1, barrier=(75.8008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59866,0.0445703,-4.03449e-05,2.09023e-08,-4.18043e-12,63457.8,21.3678], Tmin=(100,'K'), Tmax=(1443.29,'K')), NASAPolynomial(coeffs=[9.29115,0.0160375,-3.19399e-06,2.79144e-10,-8.3403e-15,61988.6,-15.9639], Tmin=(1443.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH]C(=C)C([CH2])[CH2](15717)',
    structure = SMILES('[CH]C(=C)C([CH2])[CH2]'),
    E0 = (697.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,623.023,623.023,623.023,623.023],'cm^-1')),
        HinderedRotor(inertia=(0.200176,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200177,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200177,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200177,'amu*angstrom^2'), symmetry=1, barrier=(55.1377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.85,'J/mol'), sigma=(5.90911,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.46 K, Pc=35.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21025,0.0555236,-2.63642e-05,-4.36022e-09,6.27815e-12,84020.2,26.8272], Tmin=(100,'K'), Tmax=(923.38,'K')), NASAPolynomial(coeffs=[10.335,0.0309115,-1.06122e-05,1.76033e-09,-1.15183e-13,81699.2,-19.9099], Tmin=(923.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(697.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(CO[O])[C]([CH2])[CH2](20016)',
    structure = SMILES('[CH]C(CO[O])=C([CH2])[CH2]'),
    E0 = (636.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.3724,0.0823155,-8.10295e-05,4.5409e-08,-1.06925e-11,76627,30.4615], Tmin=(100,'K'), Tmax=(1007.15,'K')), NASAPolynomial(coeffs=[11.144,0.0395355,-1.73157e-05,3.23516e-09,-2.23987e-13,74457.3,-21.5822], Tmin=(1007.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C([CH2])[CH2](9843)',
    structure = SMILES('[CH]=[C]C([CH2])[CH2]'),
    E0 = (846.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.152922,'amu*angstrom^2'), symmetry=1, barrier=(5.09302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163307,'amu*angstrom^2'), symmetry=1, barrier=(5.18453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0310369,'amu*angstrom^2'), symmetry=1, barrier=(67.385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01533,0.0427783,-3.36138e-05,9.2809e-09,2.0304e-12,101922,22.9576], Tmin=(100,'K'), Tmax=(792.169,'K')), NASAPolynomial(coeffs=[8.65034,0.0186874,-5.81895e-06,8.88198e-10,-5.43693e-14,100576,-9.37065], Tmin=(792.169,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(846.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])C([CH2])[CH2](20017)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])[CH2]'),
    E0 = (761.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,479.555],'cm^-1')),
        HinderedRotor(inertia=(0.0976643,'amu*angstrom^2'), symmetry=1, barrier=(15.9403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104629,'amu*angstrom^2'), symmetry=1, barrier=(80.6177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104631,'amu*angstrom^2'), symmetry=1, barrier=(80.6188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.494045,'amu*angstrom^2'), symmetry=1, barrier=(80.6171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.50634,'amu*angstrom^2'), symmetry=1, barrier=(80.6177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269947,0.0772992,-7.47279e-05,3.4959e-08,-5.25048e-12,91760.2,35.3693], Tmin=(100,'K'), Tmax=(923.393,'K')), NASAPolynomial(coeffs=[15.5168,0.023695,-7.86367e-06,1.27778e-09,-8.25316e-14,88414,-39.8451], Tmin=(923.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(761.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(CO[O])[C]([CH2])C(19573)',
    structure = SMILES('[CH]C(CO[O])=C([CH2])C'),
    E0 = (484.546,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371915,0.0847338,-8.89719e-05,5.76249e-08,-1.61226e-11,58403.7,30.1889], Tmin=(100,'K'), Tmax=(847.363,'K')), NASAPolynomial(coeffs=[8.54533,0.0461516,-2.06747e-05,3.8925e-09,-2.70012e-13,57018.5,-7.88952], Tmin=(847.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH]OO)C([CH2])[CH2](20018)',
    structure = SMILES('[CH]=C([CH]OO)C([CH2])[CH2]'),
    E0 = (609.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0544723,0.0816174,-7.26157e-05,2.90748e-08,-2.97654e-12,73492.3,35.7737], Tmin=(100,'K'), Tmax=(969.099,'K')), NASAPolynomial(coeffs=[17.4832,0.0247044,-8.47574e-06,1.42837e-09,-9.51428e-14,69366.6,-52.034], Tmin=(969.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Isobutyl) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](CO[O])C([CH2])=C(16267)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])CO[O]'),
    E0 = (416.859,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0020232,'amu*angstrom^2'), symmetry=1, barrier=(22.9715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0672792,'amu*angstrom^2'), symmetry=1, barrier=(87.2053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.80154,'amu*angstrom^2'), symmetry=1, barrier=(87.4049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0008,'amu*angstrom^2'), symmetry=1, barrier=(23.0104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.7978,'amu*angstrom^2'), symmetry=1, barrier=(87.319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242261,0.0815256,-7.89535e-05,4.09729e-08,-8.65428e-12,50272.9,30.093], Tmin=(100,'K'), Tmax=(1133.93,'K')), NASAPolynomial(coeffs=[14.2644,0.0320622,-1.35222e-05,2.50439e-09,-1.73095e-13,47092.9,-39.318], Tmin=(1133.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(ROOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)[CH]O[O](20019)',
    structure = SMILES('[CH2]C([CH2])C(=C)[CH]O[O]'),
    E0 = (514.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00447686,0.0778636,-7.22552e-05,3.57842e-08,-7.01388e-12,62055.6,35.9515], Tmin=(100,'K'), Tmax=(1281.59,'K')), NASAPolynomial(coeffs=[16.6816,0.0241917,-7.57225e-06,1.16722e-09,-7.19766e-14,57909.4,-48.1786], Tmin=(1281.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C([CH]O[O])C([CH2])C(19574)',
    structure = SMILES('[CH]=C([CH]O[O])C([CH2])C'),
    E0 = (556.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,911.898],'cm^-1')),
        HinderedRotor(inertia=(0.119179,'amu*angstrom^2'), symmetry=1, barrier=(2.74016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0146724,'amu*angstrom^2'), symmetry=1, barrier=(80.1073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.797168,'amu*angstrom^2'), symmetry=1, barrier=(18.3285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.797163,'amu*angstrom^2'), symmetry=1, barrier=(18.3284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11913,'amu*angstrom^2'), symmetry=1, barrier=(80.1069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0121493,0.0783964,-7.17136e-05,3.43693e-08,-6.54005e-12,67107.3,35.2083], Tmin=(100,'K'), Tmax=(1276.72,'K')), NASAPolynomial(coeffs=[17.1368,0.0247446,-8.67926e-06,1.45482e-09,-9.49526e-14,62734.6,-51.592], Tmin=(1276.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(COO)[C]([CH2])[CH2](20020)',
    structure = SMILES('[CH]C(COO)=C([CH2])[CH2]'),
    E0 = (484.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0980848,0.0861121,-7.73343e-05,3.765e-08,-7.62974e-12,58356.8,30.6806], Tmin=(100,'K'), Tmax=(1160.82,'K')), NASAPolynomial(coeffs=[13.5381,0.0397998,-1.74899e-05,3.2808e-09,-2.27794e-13,55236.5,-36.1638], Tmin=(1160.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(Allyl_P) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    E0 = (644.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (801.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (889.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (961.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1099.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1167.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (652.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (652.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (651.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (653.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (727.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (706.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (667.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (707.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (669.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (710.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (766.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (744.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (756.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (775.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (690.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1051.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (847.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1052.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (973.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (785.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (803.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (836.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (789.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (706.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (703.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['C#CCO[O](16808)', 'allyl(82)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C(C[CH][CH2])CO[O](19278)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]=C(C[O])C([CH2])[CH2](20009)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH]=C([CH][CH2])CO[O](19537)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]C([CH2])C(=[CH])CO[O](20010)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[C]=C(CO[O])C([CH2])[CH2](20011)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C(CO[O])C1CC1(20012)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH2]C1CC=C1CO[O](20013)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH2]C([CH2])C1=COOC1(19163)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C1COOCC1[CH2](18244)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.42e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['O(4)', '[CH2]C([CH2])C1=COC1(20014)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO] for rate rule [R3OO_DS;Cd_pri_rad_in;OOJ]
Euclidian distance = 2.449489742783178
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['O(4)', '[CH]=C1COCC1[CH2](18246)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.026e+11,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO_SSS;C_pri_rad_intra;OO] for rate rule [R4OO_SSS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C(CO[O])C(=C)C(19572)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH2]C(=C)C(=C)CO[O](16266)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C(COO)C([CH2])=C(20015)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]=C(CO[O])C([CH2])=C(19489)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(T)(20)', '[CH]=C(C=C)CO[O](19302)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.41156,'m^3/(mol*s)'), n=1.94471, Ea=(9.94499,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;YJ] for rate rule [Cds-OneDeH_Cds;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]CO[O](16805)', 'allyl(82)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#CCO[O](16808)', '[CH2][CH][CH2](497)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0942128,'m^3/(mol*s)'), n=2.31088, Ea=(29.6884,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cs_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]O[O](46)', 'C#CC([CH2])[CH2](9851)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0192494,'m^3/(mol*s)'), n=2.41963, Ea=(43.2169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;CsJ] for rate rule [Ct-Cs_Ct-H;CsJ-OsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O2(2)', '[CH]C(=C)C([CH2])[CH2](15717)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.96276e+06,'m^3/(mol*s)'), n=-0.119415, Ea=(1.9085,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_Sp-5R!H-4R!H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CO[O](16805)', '[CH2][CH][CH2](497)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(308.407,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH]=C(CO[O])[C]([CH2])[CH2](20016)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]O[O](46)', '[CH]=[C]C([CH2])[CH2](9843)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]=C([CH]O[O])C([CH2])[CH2](20017)'],
    products = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C(CO[O])[C]([CH2])C(19573)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C([CH]OO)C([CH2])[CH2](20018)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH2][C](CO[O])C([CH2])=C(16267)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH2]C([CH2])C(=C)[CH]O[O](20019)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C([CH]O[O])C([CH2])C(19574)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(608,'s^-1'), n=2.77, Ea=(62.3834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 324 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C(CO[O])C([CH2])[CH2](19280)'],
    products = ['[CH]=C(COO)[C]([CH2])[CH2](20020)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_SSSS;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4283',
    isomers = [
        '[CH]=C(CO[O])C([CH2])[CH2](19280)',
    ],
    reactants = [
        ('C#CCO[O](16808)', 'allyl(82)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4283',
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


species(
    label = '[CH]C=CCC(=[CH])[O](10114)',
    structure = SMILES('[CH]C=CCC(=[CH])[O]'),
    E0 = (593.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,373.022,373.029,373.061,373.225,373.278],'cm^-1')),
        HinderedRotor(inertia=(0.51927,'amu*angstrom^2'), symmetry=1, barrier=(51.2589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519496,'amu*angstrom^2'), symmetry=1, barrier=(51.2557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518521,'amu*angstrom^2'), symmetry=1, barrier=(51.2563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802791,0.0627587,-4.61968e-05,1.43259e-08,-6.16529e-13,71533.4,29.1761], Tmin=(100,'K'), Tmax=(1053.57,'K')), NASAPolynomial(coeffs=[13.7149,0.0254061,-9.63111e-06,1.7012e-09,-1.15409e-13,68164.9,-36.8647], Tmin=(1053.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2CHCHCH(4849)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6426,0.0163332,3.86245e-05,-6.71404e-08,2.83615e-11,41729.6,13.2819], Tmin=(100,'K'), Tmax=(937.716,'K')), NASAPolynomial(coeffs=[12.9703,0.0066915,-1.00084e-06,1.67635e-10,-1.71464e-14,38279.7,-43.9468], Tmin=(937.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]C=CC[C]=[CH](14989)',
    structure = SMILES('[CH]C=CC[C]=[CH]'),
    E0 = (908.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,318.762,318.776,319.035,319.411],'cm^-1')),
        HinderedRotor(inertia=(0.708805,'amu*angstrom^2'), symmetry=1, barrier=(51.0555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707921,'amu*angstrom^2'), symmetry=1, barrier=(51.0585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707738,'amu*angstrom^2'), symmetry=1, barrier=(51.0573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57614,0.0505578,-3.35032e-05,1.15347e-08,-1.66553e-12,109338,25.3154], Tmin=(100,'K'), Tmax=(1551.78,'K')), NASAPolynomial(coeffs=[10.4384,0.0277137,-1.14215e-05,2.04809e-09,-1.37184e-13,106588,-21.334], Tmin=(1551.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(908.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    label = '[C]=C([O])CC=C[CH](14990)',
    structure = SMILES('[C]C(=O)CC=C[CH]'),
    E0 = (856.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,418.05,418.055,418.055,418.056,418.056,418.058],'cm^-1')),
        HinderedRotor(inertia=(0.392749,'amu*angstrom^2'), symmetry=1, barrier=(48.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.392745,'amu*angstrom^2'), symmetry=1, barrier=(48.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.392747,'amu*angstrom^2'), symmetry=1, barrier=(48.7086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09248,0.0551949,-3.30834e-05,6.70556e-09,1.27746e-13,103080,25.0881], Tmin=(100,'K'), Tmax=(1364.31,'K')), NASAPolynomial(coeffs=[15.3678,0.0245286,-1.16668e-05,2.25059e-09,-1.57233e-13,98143.8,-52.0321], Tmin=(1364.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(856.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet) + radical(CJ3)"""),
)

species(
    label = '[CH]C=CCC1=CO1(14991)',
    structure = SMILES('[CH]C=CCC1=CO1'),
    E0 = (531.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03466,0.0528351,-9.86342e-06,-2.74554e-08,1.49324e-11,64007.5,26.1909], Tmin=(100,'K'), Tmax=(979.009,'K')), NASAPolynomial(coeffs=[15.2655,0.0233432,-8.57618e-06,1.56165e-09,-1.11025e-13,59848,-49.1752], Tmin=(979.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][CH]C[CH]C(=[CH])[O](14992)',
    structure = SMILES('[CH][CH]CC=C([CH])[O]'),
    E0 = (860.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,350,440,435,1725,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05035,0.0663716,-6.42424e-05,3.66518e-08,-8.85262e-12,103571,31.1291], Tmin=(100,'K'), Tmax=(982.421,'K')), NASAPolynomial(coeffs=[9.13535,0.033453,-1.39812e-05,2.54482e-09,-1.73341e-13,101982,-7.73302], Tmin=(982.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1[CH]CC(=[CH])O1(14993)',
    structure = SMILES('[CH]C1[CH]CC(=[CH])O1'),
    E0 = (709.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61231,0.0258625,7.81727e-05,-1.34457e-07,5.74208e-11,85388,25.4896], Tmin=(100,'K'), Tmax=(929.092,'K')), NASAPolynomial(coeffs=[22.406,0.00475943,1.78081e-06,-3.80479e-10,1.56038e-14,78571.1,-89.1909], Tmin=(929.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCJCO) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1[CH]CC([O])=C1(14994)',
    structure = SMILES('[CH]C1[CH]CC([O])=C1'),
    E0 = (548.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75836,0.0341168,2.78941e-05,-6.3281e-08,2.7361e-11,66088.8,24.9263], Tmin=(100,'K'), Tmax=(958.736,'K')), NASAPolynomial(coeffs=[14.5947,0.0181332,-5.88175e-06,1.08082e-09,-8.06902e-14,61900.7,-45.466], Tmin=(958.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=C(C)OJ) + radical(cyclopentene-4) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C1CC(=[CH])O1(14995)',
    structure = SMILES('[CH][CH]C1CC(=[CH])O1'),
    E0 = (779.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09882,0.0439892,2.55428e-05,-8.30265e-08,4.10422e-11,93908.1,24.0718], Tmin=(100,'K'), Tmax=(901.417,'K')), NASAPolynomial(coeffs=[22.0391,0.00432097,2.93636e-06,-7.68869e-10,5.23346e-14,87969.3,-86.7807], Tmin=(901.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(779.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCJCO) + radical(Cds_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C1C=C([O])C1(14996)',
    structure = SMILES('[CH][CH]C1[CH]C(=O)C1'),
    E0 = (632.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90071,0.0319683,2.96855e-05,-5.83693e-08,2.33641e-11,76142.8,26.4558], Tmin=(100,'K'), Tmax=(1003.41,'K')), NASAPolynomial(coeffs=[12.4975,0.0240306,-9.73164e-06,1.89202e-09,-1.39211e-13,72289.2,-33.3095], Tmin=(1003.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(Cs_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=CCC#C(14997)',
    structure = SMILES('[CH]C=CCC#C'),
    E0 = (592.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,262.066,262.62,262.878,263.272],'cm^-1')),
        HinderedRotor(inertia=(1.02765,'amu*angstrom^2'), symmetry=1, barrier=(50.3623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02904,'amu*angstrom^2'), symmetry=1, barrier=(50.3684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02916,'amu*angstrom^2'), symmetry=1, barrier=(50.3677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52701,0.0476137,-2.38154e-05,1.63138e-09,1.66523e-12,71390.7,21.7855], Tmin=(100,'K'), Tmax=(1137.09,'K')), NASAPolynomial(coeffs=[10.0409,0.0277972,-1.10415e-05,1.97921e-09,-1.34293e-13,68799.4,-23.2633], Tmin=(1137.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=C[CH2](8992)',
    structure = SMILES('[CH]C=C[CH2]'),
    E0 = (493.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,273.892,274.283,274.415],'cm^-1')),
        HinderedRotor(inertia=(0.949581,'amu*angstrom^2'), symmetry=1, barrier=(50.5951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.945904,'amu*angstrom^2'), symmetry=1, barrier=(50.5956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54977,0.0249465,9.89964e-06,-2.57642e-08,1.03662e-11,59460.2,14.975], Tmin=(100,'K'), Tmax=(1021.28,'K')), NASAPolynomial(coeffs=[7.69501,0.020986,-8.0647e-06,1.48637e-09,-1.04578e-13,57564.8,-14.0901], Tmin=(1021.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[CH](13797)',
    structure = SMILES('[CH]C=[CH]'),
    E0 = (623.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.17869,'amu*angstrom^2'), symmetry=1, barrier=(50.0924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.23797,0.0118291,1.62975e-05,-2.54943e-08,9.56409e-12,75051.1,10.6351], Tmin=(100,'K'), Tmax=(1001.24,'K')), NASAPolynomial(coeffs=[5.72385,0.0140485,-5.23096e-06,9.60843e-10,-6.79051e-14,73944.3,-4.4024], Tmin=(1001.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH2])[O](9653)',
    structure = SMILES('[CH]C(=C)[O]'),
    E0 = (322.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,387.097,387.511,388.12,388.534],'cm^-1')),
        HinderedRotor(inertia=(0.483854,'amu*angstrom^2'), symmetry=1, barrier=(51.3097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72438,0.0238014,-4.51958e-06,-8.49326e-09,4.40325e-12,38798.3,14.6668], Tmin=(100,'K'), Tmax=(1019.83,'K')), NASAPolynomial(coeffs=[7.15237,0.0155204,-5.70458e-06,1.01811e-09,-7.00768e-14,37422.6,-9.09941], Tmin=(1019.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][O](6861)',
    structure = SMILES('[CH][C]=O'),
    E0 = (425.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1491.4,1492.35],'cm^-1')),
        HinderedRotor(inertia=(0.0736325,'amu*angstrom^2'), symmetry=1, barrier=(5.08901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62981,0.00805288,-4.35162e-06,1.03741e-09,-8.43183e-14,51134.2,10.4631], Tmin=(100,'K'), Tmax=(1941.7,'K')), NASAPolynomial(coeffs=[6.62042,0.00292938,-1.19496e-06,2.28726e-10,-1.56225e-14,49777.3,-6.4529], Tmin=(1941.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]C=C[CH]C(=[CH])[O](14998)',
    structure = SMILES('[CH]C=CC=C([CH])[O]'),
    E0 = (686.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.715494,0.0631503,-3.47246e-05,-2.86533e-09,6.72179e-12,82691.5,27.1716], Tmin=(100,'K'), Tmax=(965.822,'K')), NASAPolynomial(coeffs=[14.4825,0.0265445,-9.57279e-06,1.6547e-09,-1.12112e-13,79080.2,-43.696], Tmin=(965.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(686.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[C]CC(=[CH])[O](14999)',
    structure = SMILES('[CH]C=[C]CC(=[CH])[O]'),
    E0 = (831.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,387.84,387.912,387.992,388.003,388.034],'cm^-1')),
        HinderedRotor(inertia=(0.486468,'amu*angstrom^2'), symmetry=1, barrier=(51.9514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486504,'amu*angstrom^2'), symmetry=1, barrier=(51.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486289,'amu*angstrom^2'), symmetry=1, barrier=(51.9515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9669,0.065029,-6.28244e-05,3.30074e-08,-7.06699e-12,100128,29.0616], Tmin=(100,'K'), Tmax=(1121.32,'K')), NASAPolynomial(coeffs=[11.9375,0.0258943,-1.04734e-05,1.88266e-09,-1.27663e-13,97667.9,-25.1213], Tmin=(1121.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CCC(=[CH])[O](15000)',
    structure = SMILES('[CH][C]=CCC(=[CH])[O]'),
    E0 = (831.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,1685,370,3120,650,792.5,1650,387.917,387.943,387.951,387.995,388.018],'cm^-1')),
        HinderedRotor(inertia=(0.486491,'amu*angstrom^2'), symmetry=1, barrier=(51.9518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486479,'amu*angstrom^2'), symmetry=1, barrier=(51.9516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486376,'amu*angstrom^2'), symmetry=1, barrier=(51.95,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.9669,0.065029,-6.28244e-05,3.30074e-08,-7.06699e-12,100128,29.0616], Tmin=(100,'K'), Tmax=(1121.32,'K')), NASAPolynomial(coeffs=[11.9375,0.0258943,-1.04734e-05,1.88266e-09,-1.27663e-13,97667.9,-25.1213], Tmin=(1121.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=C[CH]C(=[CH])O(15001)',
    structure = SMILES('[CH]C=CC=C([CH])O'),
    E0 = (548.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370072,0.0668683,-2.62527e-05,-2.13057e-08,1.53865e-11,66133.5,26.911], Tmin=(100,'K'), Tmax=(942.848,'K')), NASAPolynomial(coeffs=[18.1228,0.0233523,-7.61236e-06,1.28549e-09,-8.85463e-14,61372.5,-65.1865], Tmin=(942.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=C[CH]C(=C)[O](14183)',
    structure = SMILES('[CH]C=C[CH]C(=C)[O]'),
    E0 = (463.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,284.504,284.508,284.537,284.557,284.566,284.62],'cm^-1')),
        HinderedRotor(inertia=(0.86283,'amu*angstrom^2'), symmetry=1, barrier=(49.5536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862406,'amu*angstrom^2'), symmetry=1, barrier=(49.5567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862477,'amu*angstrom^2'), symmetry=1, barrier=(49.5573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696037,0.0627532,-3.5804e-05,-1.05197e-09,5.64588e-12,55882.5,26.3549], Tmin=(100,'K'), Tmax=(1001.71,'K')), NASAPolynomial(coeffs=[15.3617,0.0246235,-9.3035e-06,1.67377e-09,-1.16314e-13,51919.2,-49.5403], Tmin=(1001.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[C]CC(=[CH])O(15002)',
    structure = SMILES('[CH]C=[C]CC(=[CH])O'),
    E0 = (693.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,235.483,235.631,235.635,235.759],'cm^-1')),
        HinderedRotor(inertia=(1.24815,'amu*angstrom^2'), symmetry=1, barrier=(49.2656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24731,'amu*angstrom^2'), symmetry=1, barrier=(49.2625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2497,'amu*angstrom^2'), symmetry=1, barrier=(49.266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2487,'amu*angstrom^2'), symmetry=1, barrier=(49.2643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249408,0.0730979,-6.93804e-05,3.37679e-08,-6.46884e-12,83586.4,30.1379], Tmin=(100,'K'), Tmax=(1274.53,'K')), NASAPolynomial(coeffs=[17.1178,0.0201588,-7.07718e-06,1.17959e-09,-7.67356e-14,79286.4,-55.3347], Tmin=(1274.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(693.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[C]CC(=C)[O](14184)',
    structure = SMILES('[CH]C=[C]CC(=C)[O]'),
    E0 = (584.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,421.787,421.794,421.796,421.799,421.803,421.81],'cm^-1')),
        HinderedRotor(inertia=(0.41175,'amu*angstrom^2'), symmetry=1, barrier=(51.9843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411783,'amu*angstrom^2'), symmetry=1, barrier=(51.9841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411742,'amu*angstrom^2'), symmetry=1, barrier=(51.9842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894066,0.0632218,-5.20531e-05,2.30102e-08,-4.16651e-12,70414.7,28.9199], Tmin=(100,'K'), Tmax=(1308.14,'K')), NASAPolynomial(coeffs=[12.7307,0.0270281,-1.05509e-05,1.85939e-09,-1.24355e-13,67317.9,-31.3643], Tmin=(1308.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CCC(=[CH])O(15003)',
    structure = SMILES('[CH][C]=CCC(=[CH])O'),
    E0 = (693.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,235.483,235.631,235.635,235.759],'cm^-1')),
        HinderedRotor(inertia=(1.24815,'amu*angstrom^2'), symmetry=1, barrier=(49.2656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24731,'amu*angstrom^2'), symmetry=1, barrier=(49.2625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2497,'amu*angstrom^2'), symmetry=1, barrier=(49.266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2487,'amu*angstrom^2'), symmetry=1, barrier=(49.2643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249408,0.0730979,-6.93804e-05,3.37679e-08,-6.46884e-12,83586.4,30.1379], Tmin=(100,'K'), Tmax=(1274.53,'K')), NASAPolynomial(coeffs=[17.1178,0.0201588,-7.07718e-06,1.17959e-09,-7.67356e-14,79286.4,-55.3347], Tmin=(1274.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(693.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CCC(=C)[O](14185)',
    structure = SMILES('[CH][C]=CCC(=C)[O]'),
    E0 = (584.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,420.27,420.883,421.532,421.579,421.837,424.639],'cm^-1')),
        HinderedRotor(inertia=(0.413394,'amu*angstrom^2'), symmetry=1, barrier=(52.0111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411841,'amu*angstrom^2'), symmetry=1, barrier=(51.9558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409889,'amu*angstrom^2'), symmetry=1, barrier=(51.9796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894066,0.0632218,-5.20531e-05,2.30102e-08,-4.16651e-12,70414.7,28.9199], Tmin=(100,'K'), Tmax=(1308.14,'K')), NASAPolynomial(coeffs=[12.7307,0.0270281,-1.05509e-05,1.85939e-09,-1.24355e-13,67317.9,-31.3643], Tmin=(1308.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC([CH2])C(=[CH])[O](10078)',
    structure = SMILES('[CH]=CC([CH2])C(=[CH])[O]'),
    E0 = (677.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3115,3125,620,680,785,800,1600,1700,321.016],'cm^-1')),
        HinderedRotor(inertia=(0.126758,'amu*angstrom^2'), symmetry=1, barrier=(9.69332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13272,'amu*angstrom^2'), symmetry=1, barrier=(9.66036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12704,'amu*angstrom^2'), symmetry=1, barrier=(9.69409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3817.52,'J/mol'), sigma=(6.40826,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.29 K, Pc=32.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.68526,0.0740521,-9.06977e-05,5.68675e-08,-1.31771e-11,81627,28.4629], Tmin=(100,'K'), Tmax=(798.286,'K')), NASAPolynomial(coeffs=[13.0123,0.0210647,-7.63123e-06,1.27497e-09,-8.20498e-14,79379.2,-29.9831], Tmin=(798.286,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(T)(1343)',
    structure = SMILES('[CH]=[CH]'),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([597.285,1197.51,1197.81,1199.36,2685.83,3758.62],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88059,-0.000843477,2.04095e-05,-2.51946e-08,9.47585e-12,71059.3,4.72114], Tmin=(100,'K'), Tmax=(935.365,'K')), NASAPolynomial(coeffs=[4.9214,0.00358273,-9.24444e-07,1.57274e-10,-1.19541e-14,70476.2,-2.3065], Tmin=(935.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]CC(=[CH])[O](9655)',
    structure = SMILES('[CH]CC(=[CH])[O]'),
    E0 = (602.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,485.679,486.003,486.553,3552.34],'cm^-1')),
        HinderedRotor(inertia=(0.0198899,'amu*angstrom^2'), symmetry=1, barrier=(3.32905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119761,'amu*angstrom^2'), symmetry=1, barrier=(19.9624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63798,0.0445574,-4.56301e-05,2.33077e-08,-4.55552e-12,72499.5,21.5227], Tmin=(100,'K'), Tmax=(1363.28,'K')), NASAPolynomial(coeffs=[12.9679,0.00812049,-2.02506e-06,2.65656e-10,-1.49276e-14,69707.1,-35.5602], Tmin=(1363.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=C[CH]CC(=[CH])[O](15004)',
    structure = SMILES('[C]C=CCC(=[CH])[O]'),
    E0 = (892.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,215.833,215.836,215.838],'cm^-1')),
        HinderedRotor(inertia=(0.396215,'amu*angstrom^2'), symmetry=1, barrier=(13.0972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.3962,'amu*angstrom^2'), symmetry=1, barrier=(13.0972,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777681,0.0646401,-6.57004e-05,3.20522e-08,-5.67364e-12,107470,26.6097], Tmin=(100,'K'), Tmax=(1031.76,'K')), NASAPolynomial(coeffs=[15.9322,0.0144102,-5.06426e-06,8.77829e-10,-5.96716e-14,103889,-49.1739], Tmin=(1031.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(892.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(CJ3)"""),
)

species(
    label = '[CH]=C([O])CC1C=C1(15005)',
    structure = SMILES('[CH]=C([O])CC1C=C1'),
    E0 = (491.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684645,0.0650881,-6.3065e-05,3.12968e-08,-6.09143e-12,59257.3,24.4741], Tmin=(100,'K'), Tmax=(1257.25,'K')), NASAPolynomial(coeffs=[15.7428,0.0171801,-5.90735e-06,9.88647e-10,-6.48016e-14,55470.9,-51.6202], Tmin=(1257.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1CC(=[CH])O1(14933)',
    structure = SMILES('[CH]=CC1CC(=[CH])O1'),
    E0 = (506.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23161,0.0387746,4.13181e-05,-9.97702e-08,4.69062e-11,61001.9,21.8803], Tmin=(100,'K'), Tmax=(907.375,'K')), NASAPolynomial(coeffs=[22.5973,0.00313589,3.44602e-06,-8.33317e-10,5.44595e-14,54714.3,-92.4014], Tmin=(907.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC1C=C([O])C1(14887)',
    structure = SMILES('[CH]=CC1[CH]C(=O)C1'),
    E0 = (365.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9999,0.0255672,5.49781e-05,-9.04004e-08,3.63196e-11,44053.6,23.4499], Tmin=(100,'K'), Tmax=(966.885,'K')), NASAPolynomial(coeffs=[14.5969,0.019402,-6.74056e-06,1.30428e-09,-1.00002e-13,39469.8,-48.0056], Tmin=(966.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(CCJC=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC=C=C(10121)',
    structure = SMILES('[CH]=C([O])CC=C=C'),
    E0 = (399.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.670823,'amu*angstrom^2'), symmetry=1, barrier=(15.4235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671092,'amu*angstrom^2'), symmetry=1, barrier=(15.4297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71577,0.0639151,-6.04025e-05,2.91557e-08,-5.52464e-12,48194.6,27.4365], Tmin=(100,'K'), Tmax=(1289.24,'K')), NASAPolynomial(coeffs=[15.7953,0.0171291,-5.96801e-06,1.00756e-09,-6.63522e-14,44306.3,-49.1447], Tmin=(1289.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC=CC(=[CH])O(15006)',
    structure = SMILES('[CH]=CC=CC(=[CH])O'),
    E0 = (425.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.745586,0.0821772,-8.98897e-05,4.57351e-08,-8.54803e-12,51425,27.7205], Tmin=(100,'K'), Tmax=(1535.4,'K')), NASAPolynomial(coeffs=[23.4943,0.00373045,1.69317e-06,-5.18972e-10,3.9801e-14,45784.6,-93.7448], Tmin=(1535.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC=CC(=C)[O](14204)',
    structure = SMILES('[CH]=CC=CC(=C)[O]'),
    E0 = (316.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00922,'amu*angstrom^2'), symmetry=1, barrier=(23.204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00976,'amu*angstrom^2'), symmetry=1, barrier=(23.2165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139481,0.0726991,-7.37085e-05,3.61499e-08,-6.62889e-12,38255.1,26.6447], Tmin=(100,'K'), Tmax=(1536,'K')), NASAPolynomial(coeffs=[20.1187,0.00908879,-9.88883e-07,-1.31863e-11,5.83973e-15,33312.2,-75.6158], Tmin=(1536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C=CC=C(10133)',
    structure = SMILES('[CH]=C([O])C=CC=C'),
    E0 = (316.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00922,'amu*angstrom^2'), symmetry=1, barrier=(23.204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00976,'amu*angstrom^2'), symmetry=1, barrier=(23.2165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139481,0.0726991,-7.37085e-05,3.61499e-08,-6.62889e-12,38255.1,26.6447], Tmin=(100,'K'), Tmax=(1536,'K')), NASAPolynomial(coeffs=[20.1187,0.00908879,-9.88883e-07,-1.31863e-11,5.83973e-15,33312.2,-75.6158], Tmin=(1536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CCC(=[CH])O(13741)',
    structure = SMILES('[CH]=C(O)C[CH]C#C'),
    E0 = (408.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0294293,0.0757923,-8.45921e-05,4.65376e-08,-9.57453e-12,49293.5,27.6378], Tmin=(100,'K'), Tmax=(1369.34,'K')), NASAPolynomial(coeffs=[18.6609,0.0102105,-7.19137e-07,-1.54677e-10,1.96471e-14,45204.7,-64.6463], Tmin=(1369.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CCC(=C)[O](13746)',
    structure = SMILES('C#C[CH]CC(=C)[O]'),
    E0 = (299.26,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,276.078,277.484,818.478],'cm^-1')),
        HinderedRotor(inertia=(0.358614,'amu*angstrom^2'), symmetry=1, barrier=(19.7606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218786,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.304,'amu*angstrom^2'), symmetry=1, barrier=(71.3511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571635,0.066372,-6.86045e-05,3.71894e-08,-7.74866e-12,36123.8,26.5803], Tmin=(100,'K'), Tmax=(1311.25,'K')), NASAPolynomial(coeffs=[14.9484,0.0160513,-3.64584e-06,4.03532e-10,-1.8347e-14,32909.2,-44.5562], Tmin=(1311.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=CC=CC(=[CH])[O](15007)',
    structure = SMILES('[CH]=CC=CC(=[CH])[O]'),
    E0 = (563.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3115,3125,620,680,785,800,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.02891,'amu*angstrom^2'), symmetry=1, barrier=(23.6566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02851,'amu*angstrom^2'), symmetry=1, barrier=(23.6474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.017613,0.0736514,-8.2037e-05,4.35912e-08,-8.64735e-12,67964.5,26.4736], Tmin=(100,'K'), Tmax=(1399.86,'K')), NASAPolynomial(coeffs=[20.2136,0.00660669,-1.92131e-07,-1.50742e-10,1.53468e-14,63225,-74.4861], Tmin=(1399.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CCC(=[CH])[O](15008)',
    structure = SMILES('[CH]=C([O])C[CH]C#C'),
    E0 = (546.356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,2175,525,3120,650,792.5,1650,750,770,3400,2100,306.815,310.701],'cm^-1')),
        HinderedRotor(inertia=(0.00173411,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285849,'amu*angstrom^2'), symmetry=1, barrier=(19.6539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04596,'amu*angstrom^2'), symmetry=1, barrier=(70.5243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954265,0.0642867,-6.44611e-05,2.55873e-08,-2.49908e-13,65823.9,25.6249], Tmin=(100,'K'), Tmax=(821.327,'K')), NASAPolynomial(coeffs=[14.6141,0.0142632,-3.24037e-06,3.57521e-10,-1.64064e-14,63023.5,-40.9755], Tmin=(821.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Sec_Propargyl) + radical(Cds_P)"""),
)

species(
    label = 'C2H2(1342)',
    structure = SMILES('C#C'),
    E0 = (218.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,647.284,647.285,3592.16],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08745,0.00578572,8.56332e-06,-1.72824e-08,7.83594e-12,26273.8,4.4608], Tmin=(100,'K'), Tmax=(897.945,'K')), NASAPolynomial(coeffs=[5.89068,0.00209055,4.88462e-08,-5.66903e-11,4.15085e-15,25415.9,-10.7352], Tmin=(897.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.303,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=CC[CH]C(=[CH])[O](14963)',
    structure = SMILES('[CH]C([O])=CCC=[CH]'),
    E0 = (593.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,373.022,373.029,373.061,373.225,373.278],'cm^-1')),
        HinderedRotor(inertia=(0.51927,'amu*angstrom^2'), symmetry=1, barrier=(51.2589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519496,'amu*angstrom^2'), symmetry=1, barrier=(51.2557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518521,'amu*angstrom^2'), symmetry=1, barrier=(51.2563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802791,0.0627587,-4.61968e-05,1.43259e-08,-6.16529e-13,71533.4,29.1761], Tmin=(100,'K'), Tmax=(1053.57,'K')), NASAPolynomial(coeffs=[13.7149,0.0254061,-9.63111e-06,1.7012e-09,-1.15409e-13,68164.9,-36.8647], Tmin=(1053.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CCC(=[CH])[O](15009)',
    structure = SMILES('[CH]=[C]CCC(=[CH])[O]'),
    E0 = (718.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,1685,370,3115,3125,620,680,785,800,1600,1700,349.472,349.486],'cm^-1')),
        HinderedRotor(inertia=(0.106822,'amu*angstrom^2'), symmetry=1, barrier=(9.25826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106825,'amu*angstrom^2'), symmetry=1, barrier=(9.2583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10682,'amu*angstrom^2'), symmetry=1, barrier=(9.2582,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849493,0.0684658,-7.55346e-05,4.39732e-08,-1.01942e-11,86587.6,29.1902], Tmin=(100,'K'), Tmax=(1050.83,'K')), NASAPolynomial(coeffs=[13.1471,0.0216534,-8.71079e-06,1.57774e-09,-1.07779e-13,84003.1,-30.7482], Tmin=(1050.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(718.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC=[C][CH2](10125)',
    structure = SMILES('[CH]=C([O])CC=[C][CH2]'),
    E0 = (612.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,543.639,545.112],'cm^-1')),
        HinderedRotor(inertia=(0.0615582,'amu*angstrom^2'), symmetry=1, barrier=(12.9568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.563643,'amu*angstrom^2'), symmetry=1, barrier=(12.9593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858142,'amu*angstrom^2'), symmetry=1, barrier=(12.9584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738044,0.0653007,-6.40193e-05,3.22676e-08,-6.40407e-12,73778.5,29.0545], Tmin=(100,'K'), Tmax=(1229.97,'K')), NASAPolynomial(coeffs=[15.2223,0.0181961,-6.57321e-06,1.13072e-09,-7.52686e-14,70215.5,-43.8219], Tmin=(1229.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])[CH]C=C[CH2](10123)',
    structure = SMILES('[CH]C([O])=CC=C[CH2]'),
    E0 = (433.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,324.534,324.584,324.585,324.598,324.609],'cm^-1')),
        HinderedRotor(inertia=(0.679395,'amu*angstrom^2'), symmetry=1, barrier=(50.7939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679424,'amu*angstrom^2'), symmetry=1, barrier=(50.794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679391,'amu*angstrom^2'), symmetry=1, barrier=(50.794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00649,0.0541885,-1.13512e-05,-2.89508e-08,1.68159e-11,52299.3,26.5768], Tmin=(100,'K'), Tmax=(936.668,'K')), NASAPolynomial(coeffs=[15.2672,0.0223985,-7.05936e-06,1.17411e-09,-8.03082e-14,48350.8,-48.1066], Tmin=(936.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1C[CH]C([CH])O1(15010)',
    structure = SMILES('[CH][C]1C[CH]C([CH])O1'),
    E0 = (1005.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54189,0.0448756,-1.69993e-06,-3.77283e-08,2.1528e-11,121046,24.5945], Tmin=(100,'K'), Tmax=(867.095,'K')), NASAPolynomial(coeffs=[13.2491,0.0184165,-3.58269e-06,3.5871e-10,-1.70796e-14,117981,-36.1882], Tmin=(867.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(C2CsJOCs) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=C[CH]C([CH])[O](15011)',
    structure = SMILES('[CH]C=C[CH]C([CH])[O]'),
    E0 = (917.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444324,0.072001,-6.1648e-05,2.71425e-08,-4.82256e-12,110468,28.3726], Tmin=(100,'K'), Tmax=(1340.2,'K')), NASAPolynomial(coeffs=[15.5558,0.0268989,-1.11681e-05,2.03182e-09,-1.38438e-13,106418,-48.9563], Tmin=(1340.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(917.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[C]CC([CH])[O](15012)',
    structure = SMILES('[CH]C=[C]CC([CH])[O]'),
    E0 = (1038.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1685,370,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980515,0.068481,-6.411e-05,3.38454e-08,-7.55129e-12,124986,29.7266], Tmin=(100,'K'), Tmax=(1054.13,'K')), NASAPolynomial(coeffs=[10.0611,0.0340241,-1.50791e-05,2.83677e-09,-1.97248e-13,123072,-14.5606], Tmin=(1054.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1038.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Cds_S) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=[C]C[C]([CH])O(15013)',
    structure = SMILES('[CH]C=[C]C[C]([CH])O'),
    E0 = (984.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449087,0.0838297,-0.00011676,9.46744e-08,-3.11147e-11,118539,31.2446], Tmin=(100,'K'), Tmax=(826.976,'K')), NASAPolynomial(coeffs=[8.92482,0.0356007,-1.61611e-05,3.00097e-09,-2.04053e-13,117384,-6.54036], Tmin=(826.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(984.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CCC([CH])[O](15014)',
    structure = SMILES('[CH][C]=CCC([CH])[O]'),
    E0 = (1038.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1685,370,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980515,0.068481,-6.411e-05,3.38454e-08,-7.55129e-12,124986,29.7266], Tmin=(100,'K'), Tmax=(1054.13,'K')), NASAPolynomial(coeffs=[10.0611,0.0340241,-1.50791e-05,2.83677e-09,-1.97248e-13,123072,-14.5606], Tmin=(1054.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1038.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCJ2_triplet) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CC[C]([CH])O(15015)',
    structure = SMILES('[CH][C]=CC[C]([CH])O'),
    E0 = (984.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449087,0.0838297,-0.00011676,9.46744e-08,-3.11147e-11,118539,31.2446], Tmin=(100,'K'), Tmax=(826.976,'K')), NASAPolynomial(coeffs=[8.92482,0.0356007,-1.61611e-05,3.00097e-09,-2.04053e-13,117384,-6.54036], Tmin=(826.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(984.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'HCCCH(T)(14186)',
    structure = SMILES('[CH]C#C'),
    E0 = (530.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68739,0.0353212,-7.3125e-05,7.03808e-08,-2.45316e-11,63833.3,8.93421], Tmin=(100,'K'), Tmax=(908.451,'K')), NASAPolynomial(coeffs=[4.78005,0.010061,-4.92176e-06,8.86794e-10,-5.66853e-14,64115.2,2.6835], Tmin=(908.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCCH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]C(=C)O(10141)',
    structure = SMILES('[CH]C(=C)O'),
    E0 = (184.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38644,0.027409,4.34284e-06,-2.73997e-08,1.32371e-11,22240.4,14.3828], Tmin=(100,'K'), Tmax=(947.673,'K')), NASAPolynomial(coeffs=[10.6449,0.012572,-3.86221e-06,6.65059e-10,-4.73709e-14,19776.1,-29.7587], Tmin=(947.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=COC([CH])=C(15016)',
    structure = SMILES('[CH]C(=C)O[CH]C=[CH]'),
    E0 = (610.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51064,0.0665636,-5.30616e-05,2.14831e-08,-3.48127e-12,73618.1,28.7765], Tmin=(100,'K'), Tmax=(1472.57,'K')), NASAPolynomial(coeffs=[16.3799,0.0234581,-9.1537e-06,1.60527e-09,-1.06637e-13,68944.3,-53.9251], Tmin=(1472.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]1C[CH]C=CO1(15017)',
    structure = SMILES('[CH][C]1CC=C[CH]O1'),
    E0 = (583.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71185,0.0442144,-1.5726e-05,-6.34822e-09,4.69965e-12,70229.1,21.4573], Tmin=(100,'K'), Tmax=(1022.08,'K')), NASAPolynomial(coeffs=[8.74183,0.0287318,-1.06588e-05,1.86237e-09,-1.25409e-13,68163.7,-15.6856], Tmin=(1022.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(C2CsJOCs) + radical(C=CCJ(O)C) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])CC1C=[CH](14986)',
    structure = SMILES('[CH]C1([O])CC1C=[CH]'),
    E0 = (824.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00573,0.0520023,-1.07893e-05,-3.17474e-08,1.81931e-11,99262.6,25.5232], Tmin=(100,'K'), Tmax=(953.757,'K')), NASAPolynomial(coeffs=[18.1202,0.014002,-4.14675e-06,7.41109e-10,-5.58631e-14,94461.8,-64.2877], Tmin=(953.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1([O])C=C[CH]C1(15018)',
    structure = SMILES('[CH]C1([O])[CH]C=CC1'),
    E0 = (601.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20249,0.0442472,1.64393e-05,-5.91473e-08,2.71531e-11,72497.4,20.6475], Tmin=(100,'K'), Tmax=(966.473,'K')), NASAPolynomial(coeffs=[18.1377,0.0162656,-5.48708e-06,1.05873e-09,-8.18523e-14,67257.3,-70.6521], Tmin=(966.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)2OJ) + radical(C=CCJCO) + radical(CCJ2_triplet)"""),
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
    E0 = (593.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1427.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1067.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (596.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (883.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (709.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (631.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (779.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (647.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (851.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (704.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (945.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (918.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (898.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1043.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1043.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (763.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (738.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (794.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (638.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (738.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (626.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (837.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1248.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1104.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (601.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (601.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (601.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (616.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (671.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (671.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (618.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (618.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (618.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (783.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (777.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (838.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (775.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (709.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (909.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (807.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (638.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (1062.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (940.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (1127.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (1023.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (1063.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (1009.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (835.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (924.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (611.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (824.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (609.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['HCCO(2227)', 'CH2CHCHCH(4849)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH]C=CC[C]=[CH](14989)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[C]=C([O])CC=C[CH](14990)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C=CCC1=CO1(14991)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH]C[CH]C(=[CH])[O](14992)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C1[CH]CC(=[CH])O1(14993)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(115.286,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 110.5 to 115.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C1[CH]CC([O])=C1(14994)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.29e+09,'s^-1'), n=0.62, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH][CH]C1CC(=[CH])O1(14995)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(186.021,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 183.7 to 186.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH][CH]C1C=C([O])C1(14996)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O(4)', '[CH]C=CCC#C(14997)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['HCCO(2227)', '[CH]C=C[CH2](8992)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.000178913,'m^3/(mol*s)'), n=2.74787, Ea=(45.8259,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-CdHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C=[CH](13797)', '[CH]=C([CH2])[O](9653)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(212954,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C][O](6861)', '[CH]C=C[CH2](8992)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.36745e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]C=C[CH]C(=[CH])[O](14998)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH]C=[C]CC(=[CH])[O](14999)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH][C]=CCC(=[CH])[O](15000)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C=C[CH]C(=[CH])O(15001)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(587605,'s^-1'), n=2.09, Ea=(169.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C=C[CH]C(=C)[O](14183)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]C=[C]CC(=[CH])O(15002)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C=[C]CC(=C)[O](14184)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH][C]=CCC(=[CH])O(15003)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH][C]=CCC(=C)[O](14185)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC([CH2])C(=[CH])[O](10078)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C2H2(T)(1343)', '[CH]CC(=[CH])[O](9655)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[C]=C[CH]CC(=[CH])[O](15004)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C([O])CC1C=C1(15005)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=CC1CC(=[CH])O1(14933)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=CC1C=C([O])C1(14887)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C([O])CC=C=C(10121)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=CC=CC(=[CH])O(15006)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=CC=CC(=C)[O](14204)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C([O])C=CC=C(10133)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C=CCC(=[CH])O(13741)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C=CCC(=C)[O](13746)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH]=CC=CC(=[CH])[O](15007)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH]=C=CCC(=[CH])[O](15008)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1476.16,'m^3/(mol*s)'), n=1.53926, Ea=(19.2512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cdd_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C2H2(1342)', '[CH]CC(=[CH])[O](9655)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3772.39,'m^3/(mol*s)'), n=1.05332, Ea=(17.9296,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-H;YJ] for rate rule [Ct-H_Ct-H;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C][O](6861)', 'CH2CHCHCH(4849)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.0105759,'m^3/(mol*s)'), n=2.40084, Ea=(4.47938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OneDeH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=CC[CH]C(=[CH])[O](14963)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CCC(=[CH])[O](15009)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(9.57388e+10,'s^-1'), n=0.905639, Ea=(190.563,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C([O])CC=[C][CH2](10125)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C([O])[CH]C=C[CH2](10123)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4HJ_2;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH][C]1C[CH]C([CH])O1(15010)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]C=C[CH]C([CH])[O](15011)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C=[C]CC([CH])[O](15012)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]C=[C]C[C]([CH])O(15013)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH][C]=CCC([CH])[O](15014)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH][C]=CC[C]([CH])O(15015)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['HCCCH(T)(14186)', '[CH]C(=C)O(10141)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3.59298e+09,'s^-1'), n=0.614601, Ea=(241.463,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.07280803979615884, var=3.9519504359180297, Tref=1000.0, N=5, correlation='Root_1R!H->C_2R!H->C_5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_5R!H->O
    Total Standard Deviation in ln(k): 4.168250509264941
Exact match found for rate rule [Root_1R!H->C_2R!H->C_5R!H->O]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]C=COC([CH])=C(15016)'],
    products = ['[CH]C=CCC(=[CH])[O](10114)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH][C]1C[CH]C=CO1(15017)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4.59219e+10,'s^-1'), n=0.253963, Ea=(17.5802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R6_linear;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C1([O])CC1C=[CH](14986)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(230.564,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHCd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 229.0 to 230.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]C1([O])C=C[CH]C1(15018)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.92551e+11,'s^-1'), n=0.201102, Ea=(15.495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R6;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = 'PDepNetwork #3350',
    isomers = [
        '[CH]C=CCC(=[CH])[O](10114)',
    ],
    reactants = [
        ('HCCO(2227)', 'CH2CHCHCH(4849)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3350',
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


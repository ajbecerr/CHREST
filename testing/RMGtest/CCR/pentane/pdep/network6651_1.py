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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.379689,0.109135,-0.000186942,1.6975e-07,-5.90408e-11,80777.5,33.8989], Tmin=(100,'K'), Tmax=(861.112,'K')), NASAPolynomial(coeffs=[8.68495,0.0400147,-1.94837e-05,3.67521e-09,-2.49097e-13,80217.9,-2.66198], Tmin=(861.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S)"""),
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
    label = '[CH2]C([O])(C=O)C[C]=C(24963)',
    structure = SMILES('[CH2]C([O])(C=O)C[C]=C'),
    E0 = (364.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180,1600,1804.55,2711.8,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159362,'amu*angstrom^2'), symmetry=1, barrier=(3.66404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159362,'amu*angstrom^2'), symmetry=1, barrier=(3.66404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159362,'amu*angstrom^2'), symmetry=1, barrier=(3.66404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159362,'amu*angstrom^2'), symmetry=1, barrier=(3.66404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0553746,0.0945316,-0.000136482,1.07466e-07,-3.3692e-11,43985.6,33.4616], Tmin=(100,'K'), Tmax=(829.114,'K')), NASAPolynomial(coeffs=[12.2312,0.0309696,-1.37334e-05,2.53197e-09,-1.71414e-13,42095.6,-22.6229], Tmin=(829.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)(C=O)OJ) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])C([O])C[C]=C(25957)',
    structure = SMILES('[CH2][C]([O])C([O])C[C]=C'),
    E0 = (677.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0632773,0.093146,-0.000125745,9.0191e-08,-2.48364e-11,81651.1,33.5831], Tmin=(100,'K'), Tmax=(693.501,'K')), NASAPolynomial(coeffs=[12.0899,0.0329031,-1.51785e-05,2.87461e-09,-1.98809e-13,79763.6,-21.6184], Tmin=(693.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC[C]([O])[CH][O](25975)',
    structure = SMILES('C=[C]CC[C]([O])[CH][O]'),
    E0 = (654.453,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4355.37,'J/mol'), sigma=(7.29009,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.30 K, Pc=25.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0459904,0.102786,-0.000177704,1.66511e-07,-5.95273e-11,78845,34.2652], Tmin=(100,'K'), Tmax=(859.031,'K')), NASAPolynomial(coeffs=[6.14375,0.0436093,-2.13677e-05,4.04847e-09,-2.75335e-13,78901.5,11.8627], Tmin=(859.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CCsJOH) + radical(Cds_S)"""),
)

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
    label = '[CH2][C]([CH][O])C[C]=C(27761)',
    structure = SMILES('[CH2][C]([CH][O])C[C]=C'),
    E0 = (768.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180,818.174,829.105,830.649],'cm^-1')),
        HinderedRotor(inertia=(0.148663,'amu*angstrom^2'), symmetry=1, barrier=(3.41805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00681977,'amu*angstrom^2'), symmetry=1, barrier=(3.2776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147758,'amu*angstrom^2'), symmetry=1, barrier=(3.39726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143705,'amu*angstrom^2'), symmetry=1, barrier=(3.30406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00968,0.0725356,-0.000101231,8.77742e-08,-3.0668e-11,92536.3,30.5849], Tmin=(100,'K'), Tmax=(843.728,'K')), NASAPolynomial(coeffs=[5.51945,0.0382102,-1.71917e-05,3.18674e-09,-2.16229e-13,92236.1,12.3247], Tmin=(843.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(768.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CCsJOH) + radical(Isobutyl) + radical(Cds_S)"""),
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
    label = '[CH2][C]CC(=C)[O](2630)',
    structure = SMILES('[CH2][C]CC(=C)[O]'),
    E0 = (536.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180,180.005,813.639,815.658],'cm^-1')),
        HinderedRotor(inertia=(0.637143,'amu*angstrom^2'), symmetry=1, barrier=(14.6492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684077,'amu*angstrom^2'), symmetry=1, barrier=(1.57283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311223,'amu*angstrom^2'), symmetry=1, barrier=(14.6464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11353,0.0564151,-5.27537e-05,2.54942e-08,-4.85736e-12,64624.6,25.9357], Tmin=(100,'K'), Tmax=(1279.99,'K')), NASAPolynomial(coeffs=[13.987,0.0161851,-5.60879e-06,9.3933e-10,-6.14419e-14,61329,-39.3494], Tmin=(1279.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][C]CC([O])=C[O](12327)',
    structure = SMILES('[CH2][C]CC([O])=C[O]'),
    E0 = (469.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,360.28,360.368,360.528,360.529,360.772],'cm^-1')),
        HinderedRotor(inertia=(0.00129779,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236601,'amu*angstrom^2'), symmetry=1, barrier=(21.8177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236389,'amu*angstrom^2'), symmetry=1, barrier=(21.8137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.025734,0.0711809,-7.44581e-05,3.70155e-08,-6.8608e-12,56576.6,31.2999], Tmin=(100,'K'), Tmax=(1513.02,'K')), NASAPolynomial(coeffs=[20.4528,0.00678045,-4.3938e-07,-8.07887e-11,9.3185e-15,51554.2,-72.0965], Tmin=(1513.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'H2CC(T)(1341)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.3669e-08,8.45456e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.859e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([CH2])([O])[CH][O](26275)',
    structure = SMILES('[CH2]C([CH2])([O])[CH][O]'),
    E0 = (564.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,180,180,180,280.802,1488.31,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.13832,'amu*angstrom^2'), symmetry=1, barrier=(3.18024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13832,'amu*angstrom^2'), symmetry=1, barrier=(3.18024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13832,'amu*angstrom^2'), symmetry=1, barrier=(3.18024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440529,0.0919539,-0.00017559,1.64977e-07,-5.75825e-11,67992.7,24.6642], Tmin=(100,'K'), Tmax=(884.555,'K')), NASAPolynomial(coeffs=[7.32782,0.0288663,-1.44403e-05,2.69762e-09,-1.79449e-13,68023.9,-0.654324], Tmin=(884.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C([O])([C][O])C[C]=C(27762)',
    structure = SMILES('[CH2]C([O])([C][O])C[C]=C'),
    E0 = (951.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.339815,0.106539,-0.000181581,1.61622e-07,-5.54892e-11,114537,32.9385], Tmin=(100,'K'), Tmax=(848.7,'K')), NASAPolynomial(coeffs=[10.4215,0.0347164,-1.73448e-05,3.31386e-09,-2.2657e-13,113470,-12.7358], Tmin=(848.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(951.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(Cds_S) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C([O])([CH][O])C[C]=C(27763)',
    structure = SMILES('[CH]C([O])([CH][O])C[C]=C'),
    E0 = (906.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.138801,0.101814,-0.000170447,1.51869e-07,-5.22423e-11,109179,34.1998], Tmin=(100,'K'), Tmax=(854.89,'K')), NASAPolynomial(coeffs=[9.32429,0.0361063,-1.75548e-05,3.3171e-09,-2.25459e-13,108344,-5.39038], Tmin=(854.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(906.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1(C[C]=C)OC1[O](24189)',
    structure = SMILES('[CH2]C1(C[C]=C)OC1[O]'),
    E0 = (408.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4095.79,'J/mol'), sigma=(6.91196,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=639.75 K, Pc=28.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0402812,0.0882552,-0.000114345,8.13439e-08,-2.25853e-11,49315.1,29.4963], Tmin=(100,'K'), Tmax=(988.31,'K')), NASAPolynomial(coeffs=[12.8026,0.0277934,-9.20962e-06,1.40655e-09,-8.33133e-14,47222.7,-29.7479], Tmin=(988.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC1([CH][O])CO1(27764)',
    structure = SMILES('C=[C]CC1([CH][O])CO1'),
    E0 = (419.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0977711,0.0887364,-0.000119228,8.87661e-08,-2.56475e-11,50537.2,29.3245], Tmin=(100,'K'), Tmax=(980.441,'K')), NASAPolynomial(coeffs=[11.4553,0.0298981,-1.00837e-05,1.54598e-09,-9.12474e-14,48911,-22.1805], Tmin=(980.441,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC1([O])CC1[O](27765)',
    structure = SMILES('C=[C]CC1([O])CC1[O]'),
    E0 = (432.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442657,0.0673224,-4.27537e-05,2.91151e-09,4.77443e-12,52180.4,29.0228], Tmin=(100,'K'), Tmax=(1014.75,'K')), NASAPolynomial(coeffs=[17.3947,0.0223097,-8.45564e-06,1.55923e-09,-1.10653e-13,47617.1,-58.5422], Tmin=(1014.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1([CH][O])CC(=C)O1(27750)',
    structure = SMILES('[CH2]C1([CH][O])CC(=C)O1'),
    E0 = (323.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.146829,0.0851246,-8.33506e-05,3.24972e-08,-1.05295e-12,39034.5,24.5811], Tmin=(100,'K'), Tmax=(874.602,'K')), NASAPolynomial(coeffs=[19.44,0.0172548,-4.18521e-06,5.36292e-10,-3.02407e-14,34778,-72.037], Tmin=(874.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])CC(=C)C1[O](27766)',
    structure = SMILES('[CH2]C1([O])CC(=C)C1[O]'),
    E0 = (400.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296173,0.0709487,-5.28697e-05,1.54594e-08,-4.53006e-13,48283.1,27.5235], Tmin=(100,'K'), Tmax=(1112.19,'K')), NASAPolynomial(coeffs=[17.5371,0.0239999,-9.85914e-06,1.85149e-09,-1.30556e-13,43516.7,-61.674], Tmin=(1112.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'C=C1CC([O])([CH][O])C1(27642)',
    structure = SMILES('C=C1CC([O])([CH][O])C1'),
    E0 = (376.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15639,0.0677523,-6.03047e-05,3.12609e-08,-7.13944e-12,45425.2,25.3395], Tmin=(100,'K'), Tmax=(1003.31,'K')), NASAPolynomial(coeffs=[8.12833,0.0399561,-1.87474e-05,3.64705e-09,-2.58648e-13,44026.2,-8.31903], Tmin=(1003.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(O)([CH][O])C=C=C(27767)',
    structure = SMILES('[CH2]C(O)([CH][O])C=C=C'),
    E0 = (362.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.155765,0.0988021,-0.000150971,1.25551e-07,-4.1083e-11,43687,33.639], Tmin=(100,'K'), Tmax=(850.843,'K')), NASAPolynomial(coeffs=[11.1348,0.0338056,-1.53762e-05,2.84545e-09,-1.92081e-13,42197.1,-16.4729], Tmin=(850.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([O])(C=C=C)C[O](27768)',
    structure = SMILES('[CH2]C([O])(C=C=C)C[O]'),
    E0 = (410.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51821,0.0835142,-0.000115873,9.84456e-08,-3.44535e-11,49533,32.8844], Tmin=(100,'K'), Tmax=(791.903,'K')), NASAPolynomial(coeffs=[7.20891,0.0411864,-1.95349e-05,3.73772e-09,-2.59412e-13,48740.9,3.85588], Tmin=(791.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'C=C=CC(C)([O])[CH][O](27769)',
    structure = SMILES('C=C=CC(C)([O])[CH][O]'),
    E0 = (377.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400104,0.0859403,-0.000120224,1.00706e-07,-3.45609e-11,45551,31.9828], Tmin=(100,'K'), Tmax=(796.179,'K')), NASAPolynomial(coeffs=[8.0619,0.0397523,-1.87083e-05,3.56388e-09,-2.46618e-13,44574.8,-1.70284], Tmin=(796.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])([CH][C]=C)C=O(25060)',
    structure = SMILES('[CH2]C([O])([CH][C]=C)C=O'),
    E0 = (481.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,180,180,180,297.629,1440.72,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.139315,'amu*angstrom^2'), symmetry=1, barrier=(3.20314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139315,'amu*angstrom^2'), symmetry=1, barrier=(3.20314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139315,'amu*angstrom^2'), symmetry=1, barrier=(3.20314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139315,'amu*angstrom^2'), symmetry=1, barrier=(3.20314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.229211,0.0979098,-0.000141946,1.07917e-07,-3.24232e-11,58054.2,31.288], Tmin=(100,'K'), Tmax=(818.381,'K')), NASAPolynomial(coeffs=[14.2485,0.0271516,-1.22629e-05,2.28146e-09,-1.55523e-13,55684.4,-35.6576], Tmin=(818.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)(C=O)OJ) + radical(C=CCJCO) + radical(CJC(C)(C=O)O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC([CH2])([O])C=O(25127)',
    structure = SMILES('[CH]=[C]CC([CH2])([O])C=O'),
    E0 = (611.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650,180,180,180,1234.84,1600,2151.98,3200],'cm^-1')),
        HinderedRotor(inertia=(0.145633,'amu*angstrom^2'), symmetry=1, barrier=(3.34839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145633,'amu*angstrom^2'), symmetry=1, barrier=(3.34839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145633,'amu*angstrom^2'), symmetry=1, barrier=(3.34839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145633,'amu*angstrom^2'), symmetry=1, barrier=(3.34839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.170422,0.0985632,-0.00015501,1.2737e-07,-4.07083e-11,73707.1,34.2761], Tmin=(100,'K'), Tmax=(865.482,'K')), NASAPolynomial(coeffs=[12.6541,0.0278177,-1.25118e-05,2.28826e-09,-1.52779e-13,71917,-23.2596], Tmin=(865.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(611.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)(C=O)OJ) + radical(CJC(C)(C=O)O) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]CC(=C)[O](15511)',
    structure = SMILES('C=[C]CC(=C)[O]'),
    E0 = (249.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,400.2,400.397,400.974],'cm^-1')),
        HinderedRotor(inertia=(0.0759882,'amu*angstrom^2'), symmetry=1, barrier=(8.66216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0765012,'amu*angstrom^2'), symmetry=1, barrier=(8.68948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.78,'J/mol'), sigma=(6.06207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.87 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54215,0.0477417,-3.21906e-05,6.04401e-09,1.69285e-12,30142.6,23.6773], Tmin=(100,'K'), Tmax=(1008.06,'K')), NASAPolynomial(coeffs=[12.024,0.0183644,-6.65273e-06,1.17508e-09,-8.06668e-14,27408.7,-30.054], Tmin=(1008.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC([O])=C[O](15579)',
    structure = SMILES('C=[C]CC([O])=C[O]'),
    E0 = (182.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,317.264,317.319,317.333,317.341],'cm^-1')),
        HinderedRotor(inertia=(0.225433,'amu*angstrom^2'), symmetry=1, barrier=(16.1111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225583,'amu*angstrom^2'), symmetry=1, barrier=(16.1105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4284.46,'J/mol'), sigma=(6.81655,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.22 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.950635,0.0561283,-3.18828e-05,-1.06961e-08,1.16732e-11,22070.5,27.0702], Tmin=(100,'K'), Tmax=(928.985,'K')), NASAPolynomial(coeffs=[18.2946,0.0094926,-1.86303e-06,2.55946e-10,-1.89319e-14,17637.9,-61.84], Tmin=(928.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(=C)[CH][O](22312)',
    structure = SMILES('[CH2]C(=C[O])C[C]=C'),
    E0 = (371.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,353.512,353.512,353.514],'cm^-1')),
        HinderedRotor(inertia=(0.250934,'amu*angstrom^2'), symmetry=1, barrier=(22.2534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250934,'amu*angstrom^2'), symmetry=1, barrier=(22.2535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250934,'amu*angstrom^2'), symmetry=1, barrier=(22.2534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.748909,0.0587907,-2.34523e-05,-1.93985e-08,1.37988e-11,44834.4,27.2728], Tmin=(100,'K'), Tmax=(959.675,'K')), NASAPolynomial(coeffs=[18.0392,0.0169981,-5.44976e-06,9.6808e-10,-7.02382e-14,40121.7,-62.6949], Tmin=(959.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH2]C([O])([CH][O])[CH][C]=C(25056)',
    structure = SMILES('[CH2]C([O])([CH][O])[CH][C]=C'),
    E0 = (787.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.592707,0.113023,-0.000194434,1.73225e-07,-5.92604e-11,94847.8,31.8628], Tmin=(100,'K'), Tmax=(856.735,'K')), NASAPolynomial(coeffs=[10.7935,0.0360323,-1.79144e-05,3.40066e-09,-2.31166e-13,93771.3,-16.2045], Tmin=(856.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(787.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC([CH2])([O])[CH][O](25120)',
    structure = SMILES('[CH]=[C]CC([CH2])([O])[CH][O]'),
    E0 = (917.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,180,180,365.963,799.219,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154176,'amu*angstrom^2'), symmetry=1, barrier=(3.54481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154176,'amu*angstrom^2'), symmetry=1, barrier=(3.54481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154176,'amu*angstrom^2'), symmetry=1, barrier=(3.54481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154176,'amu*angstrom^2'), symmetry=1, barrier=(3.54481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.468127,0.112827,-0.000204162,1.87777e-07,-6.51719e-11,110498,34.6195], Tmin=(100,'K'), Tmax=(872.803,'K')), NASAPolynomial(coeffs=[9.03079,0.037,-1.83438e-05,3.45128e-09,-2.32131e-13,110070,-2.86837], Tmin=(872.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(917.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=CC([CH2])([O])[CH][O](27195)',
    structure = SMILES('[CH2]C([O])([CH][O])[CH]C=C'),
    E0 = (549.488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.521328,0.10795,-0.000167909,1.41635e-07,-4.73045e-11,66242.8,31.1345], Tmin=(100,'K'), Tmax=(821.64,'K')), NASAPolynomial(coeffs=[11.799,0.0368558,-1.7828e-05,3.39763e-09,-2.33703e-13,64593.4,-23.6007], Tmin=(821.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH]=CCC([CH2])([O])[CH][O](27538)',
    structure = SMILES('[CH]=CCC([CH2])([O])[CH][O]'),
    E0 = (679.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,364.295,796.24,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153963,'amu*angstrom^2'), symmetry=1, barrier=(3.53991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.41755,0.108013,-0.000178594,1.57479e-07,-5.37744e-11,81893.8,33.9648], Tmin=(100,'K'), Tmax=(853.569,'K')), NASAPolynomial(coeffs=[10.1224,0.0376712,-1.8167e-05,3.42645e-09,-2.32832e-13,80857.6,-10.7451], Tmin=(853.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O)([CH][O])[CH][C]=C(27770)',
    structure = SMILES('[CH2]C(O)([CH][O])[CH][C]=C'),
    E0 = (558.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.812971,0.117009,-0.000196548,1.70991e-07,-5.72562e-11,67324.3,32.5928], Tmin=(100,'K'), Tmax=(863.342,'K')), NASAPolynomial(coeffs=[12.0995,0.0354332,-1.70271e-05,3.18623e-09,-2.14794e-13,65905.3,-23.1107], Tmin=(863.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])([CH][C]=C)C[O](27771)',
    structure = SMILES('[CH2]C([O])([CH][C]=C)C[O]'),
    E0 = (607.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.321328,0.104261,-0.000163201,1.41076e-07,-4.82672e-11,73155.8,31.34], Tmin=(100,'K'), Tmax=(824.389,'K')), NASAPolynomial(coeffs=[10.1383,0.0392963,-1.91336e-05,3.65689e-09,-2.51892e-13,71914.3,-14.1722], Tmin=(824.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CJC(C)2O) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C(C)([O])[CH][O](27772)',
    structure = SMILES('C=[C][CH]C(C)([O])[CH][O]'),
    E0 = (575.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.260027,0.102175,-0.000156778,1.33629e-07,-4.53094e-11,69339.5,30.8311], Tmin=(100,'K'), Tmax=(820.021,'K')), NASAPolynomial(coeffs=[10.3399,0.0385463,-1.85771e-05,3.54164e-09,-2.43898e-13,68002,-15.7598], Tmin=(820.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])([CH]O)[CH][C]=C(27773)',
    structure = SMILES('[CH2]C([O])([CH]O)[CH][C]=C'),
    E0 = (561.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.921044,0.115339,-0.000180471,1.45708e-07,-4.59816e-11,67718.3,32.1828], Tmin=(100,'K'), Tmax=(843.362,'K')), NASAPolynomial(coeffs=[15.4213,0.0297276,-1.37962e-05,2.56532e-09,-1.7338e-13,65249.9,-42.1677], Tmin=(843.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC([CH2])(O)[CH][O](27774)',
    structure = SMILES('[CH]=[C]CC([CH2])(O)[CH][O]'),
    E0 = (688.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,180,180,180,1598.77,1600,2923.92,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149569,'amu*angstrom^2'), symmetry=1, barrier=(3.43889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149569,'amu*angstrom^2'), symmetry=1, barrier=(3.43889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149569,'amu*angstrom^2'), symmetry=1, barrier=(3.43889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149569,'amu*angstrom^2'), symmetry=1, barrier=(3.43889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149569,'amu*angstrom^2'), symmetry=1, barrier=(3.43889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.688628,0.116819,-0.000206311,1.85621e-07,-6.3221e-11,82974.4,35.3501], Tmin=(100,'K'), Tmax=(879.127,'K')), NASAPolynomial(coeffs=[10.3279,0.0364163,-1.74654e-05,3.23896e-09,-2.15935e-13,82207.5,-9.72475], Tmin=(879.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC([CH2])([O])C[O](27775)',
    structure = SMILES('[CH]=[C]CC([CH2])([O])C[O]'),
    E0 = (737.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,180,180,419.685,744.245,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156102,'amu*angstrom^2'), symmetry=1, barrier=(3.58909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156102,'amu*angstrom^2'), symmetry=1, barrier=(3.58909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156102,'amu*angstrom^2'), symmetry=1, barrier=(3.58909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156102,'amu*angstrom^2'), symmetry=1, barrier=(3.58909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.212697,0.104261,-0.000173635,1.56545e-07,-5.45519e-11,88806.6,34.1533], Tmin=(100,'K'), Tmax=(852.482,'K')), NASAPolynomial(coeffs=[8.45143,0.0401302,-1.94837e-05,3.68841e-09,-2.51249e-13,88182.5,-1.25989], Tmin=(852.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(737.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(C)([O])[CH][O](27776)',
    structure = SMILES('[CH]=[C]CC(C)([O])[CH][O]'),
    E0 = (705.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,180,180,180,460.53,701.408,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157513,'amu*angstrom^2'), symmetry=1, barrier=(3.62154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157513,'amu*angstrom^2'), symmetry=1, barrier=(3.62154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157513,'amu*angstrom^2'), symmetry=1, barrier=(3.62154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157513,'amu*angstrom^2'), symmetry=1, barrier=(3.62154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.155548,0.102228,-0.000167418,1.49394e-07,-5.17332e-11,84990.5,33.659], Tmin=(100,'K'), Tmax=(851.541,'K')), NASAPolynomial(coeffs=[8.66531,0.0393582,-1.89141e-05,3.56997e-09,-2.42986e-13,84265.3,-2.91578], Tmin=(851.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(705.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC([CH2])([O])[CH]O(27777)',
    structure = SMILES('[CH]=[C]CC([CH2])([O])[CH]O'),
    E0 = (691.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,180,180,180,1600,1604.19,2908.4,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149575,'amu*angstrom^2'), symmetry=1, barrier=(3.43901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149575,'amu*angstrom^2'), symmetry=1, barrier=(3.43901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149575,'amu*angstrom^2'), symmetry=1, barrier=(3.43901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149575,'amu*angstrom^2'), symmetry=1, barrier=(3.43901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149575,'amu*angstrom^2'), symmetry=1, barrier=(3.43901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.813821,0.115365,-0.000191059,1.61505e-07,-5.24867e-11,83369.1,35.0007], Tmin=(100,'K'), Tmax=(873.88,'K')), NASAPolynomial(coeffs=[13.7038,0.0306148,-1.41776e-05,2.60434e-09,-1.73366e-13,81530.5,-29.0836], Tmin=(873.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (670.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (670.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (841.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (827.731,'kJ/mol'),
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
    E0 = (1287.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1022.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (906.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1220.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1162.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1118.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (673.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (675.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (675.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (678.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (678.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (678.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (748.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (748.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (748.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (701.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (836.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (703.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (670.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (670.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (670.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (764.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (967.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (999.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1129.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (828.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (785.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (745.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (786.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (788.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (793.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (721.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (770.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (738.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (843.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=C([O])[CH][O](2850)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C([O])(C=O)C[C]=C(24963)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][C]([O])C([O])C[C]=C(25957)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=[C]CC[C]([O])[CH][O](25975)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
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
    reactants = ['O(4)', '[CH2][C]([CH][O])C[C]=C(27761)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][O](1548)', '[CH2][C]CC(=C)[O](2630)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH2][C]CC([O])=C[O](12327)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2CC(T)(1341)', '[CH2]C([CH2])([O])[CH][O](26275)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2]C([O])([C][O])C[C]=C(27762)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]C([O])([CH][O])C[C]=C(27763)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C1(C[C]=C)OC1[O](24189)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=[C]CC1([CH][O])CO1(27764)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=[C]CC1([O])CC1[O](27765)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C1([CH][O])CC(=C)O1(27750)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C1([O])CC(=C)C1[O](27766)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=C1CC([O])([CH][O])C1(27642)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C(O)([CH][O])C=C=C(27767)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C([O])(C=C=C)C[O](27768)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=C=CC(C)([O])[CH][O](27769)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2]C([O])([CH][C]=C)C=O(25060)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH]=[C]CC([CH2])([O])C=O(25127)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH][O](1548)', 'C=[C]CC(=C)[O](15511)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.75465,'m^3/(mol*s)'), n=1.70829, Ea=(22.0129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [CO_O;YJ] for rate rule [CO_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(20)', 'C=[C]CC([O])=C[O](15579)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(105.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;Y_1centerbirad] for rate rule [CO_O;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.8 to 105.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)', 'C=[C]CC(=C)[CH][O](22312)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(55.7026,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 55.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C([O])[CH][O](2850)', '[CH2][C]=C(459)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4679.9,'m^3/(mol*s)'), n=0.573452, Ea=(254.146,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 253.5 to 254.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]([O])[CH][O](26240)', 'allene(458)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00859497,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]([O])[CH][O](26240)', '[CH2][C]=C(459)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2]C([O])([CH][O])[CH][C]=C(25056)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH]=[C]CC([CH2])([O])[CH][O](25120)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C=CC([CH2])([O])[CH][O](27195)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CCC([CH2])([O])[CH][O](27538)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C(O)([CH][O])[CH][C]=C(27770)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C([O])([CH][C]=C)C[O](27771)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['C=[C][CH]C(C)([O])[CH][O](27772)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    products = ['[CH2]C([O])([CH]O)[CH][C]=C(27773)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.04154e+06,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]CC([CH2])(O)[CH][O](27774)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]CC([CH2])([O])C[O](27775)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]CC(C)([O])[CH][O](27776)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC([CH2])([O])[CH]O(27777)'],
    products = ['[CH2]C([O])([CH][O])C[C]=C(24956)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #6651',
    isomers = [
        '[CH2]C([O])([CH][O])C[C]=C(24956)',
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
    label = 'PDepNetwork #6651',
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


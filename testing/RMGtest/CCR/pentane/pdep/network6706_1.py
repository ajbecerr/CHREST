species(
    label = '[CH2][C]=CC([O])[C]([CH2])[O](25999)',
    structure = SMILES('[CH2][C]([O])C([O])[CH][C]=C'),
    E0 = (794.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.210029,0.0979404,-0.000137553,1.01495e-07,-2.97678e-11,95723.9,31.7523], Tmin=(100,'K'), Tmax=(835.745,'K')), NASAPolynomial(coeffs=[14.221,0.028867,-1.35714e-05,2.58989e-09,-1.79954e-13,93312,-35.2787], Tmin=(835.745,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(794.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C=CCJCO) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C3H3(5450)',
    structure = SMILES('[CH]=C=C'),
    E0 = (338.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,180,180,2603.58],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][C]([O])C([O])C=C=C(25984)',
    structure = SMILES('[CH2][C]([O])C([O])C=C=C'),
    E0 = (603.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,540,610,2055,415.089,415.09,415.091,415.099],'cm^-1')),
        HinderedRotor(inertia=(0.110823,'amu*angstrom^2'), symmetry=1, barrier=(13.5495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110824,'amu*angstrom^2'), symmetry=1, barrier=(13.5495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110811,'amu*angstrom^2'), symmetry=1, barrier=(13.5496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.250046,0.0851742,-0.00010329,6.52878e-08,-1.64659e-11,72769.8,31.9522], Tmin=(100,'K'), Tmax=(965.633,'K')), NASAPolynomial(coeffs=[14.4254,0.0264556,-1.20785e-05,2.31675e-09,-1.6302e-13,70032.1,-35.94], Tmin=(965.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
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
    label = '[CH2][C]([O])[CH][CH][C]=C(28221)',
    structure = SMILES('[CH2][C][CH]C=C([CH2])[O]'),
    E0 = (800.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,1051.7,1051.76],'cm^-1')),
        HinderedRotor(inertia=(0.13137,'amu*angstrom^2'), symmetry=1, barrier=(18.7546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111422,'amu*angstrom^2'), symmetry=1, barrier=(87.4532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815765,'amu*angstrom^2'), symmetry=1, barrier=(18.756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00781212,'amu*angstrom^2'), symmetry=1, barrier=(87.4676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633182,0.0673047,-6.77673e-05,3.51898e-08,-7.17149e-12,96394,29.7393], Tmin=(100,'K'), Tmax=(1202.38,'K')), NASAPolynomial(coeffs=[15.4837,0.0179007,-6.1341e-06,1.01672e-09,-6.61349e-14,92822.8,-44.6428], Tmin=(1202.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=C(O)CJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = 'C=[C][CH][CH][O](17610)',
    structure = SMILES('[CH2][C]=C[CH][O]'),
    E0 = (549.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,180,1441.42],'cm^-1')),
        HinderedRotor(inertia=(0.0637222,'amu*angstrom^2'), symmetry=1, barrier=(39.9635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270909,'amu*angstrom^2'), symmetry=1, barrier=(39.9484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50002,0.0283298,-1.09144e-05,-2.0554e-09,1.67008e-12,66123.7,20.3322], Tmin=(100,'K'), Tmax=(1235.14,'K')), NASAPolynomial(coeffs=[8.2763,0.0176346,-7.655e-06,1.43668e-09,-9.96449e-14,64085.8,-11.2287], Tmin=(1235.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH2](16918)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,260.76,263.05,263.653],'cm^-1')),
        HinderedRotor(inertia=(1.04394,'amu*angstrom^2'), symmetry=1, barrier=(50.8215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2703,0.0129731,8.09839e-06,-1.37802e-08,4.71189e-12,73935.1,11.2848], Tmin=(100,'K'), Tmax=(1104.68,'K')), NASAPolynomial(coeffs=[4.54119,0.0159815,-6.32008e-06,1.15742e-09,-7.99393e-14,73190,2.92522], Tmin=(1104.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH2][C]([O])C([O])[C][C]=C(28222)',
    structure = SMILES('[CH2][C]=[C]C([O])[C]([CH2])[O]'),
    E0 = (1054.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,362.035,363.239,365.077,2337.35],'cm^-1')),
        HinderedRotor(inertia=(0.156222,'amu*angstrom^2'), symmetry=1, barrier=(14.6687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155229,'amu*angstrom^2'), symmetry=1, barrier=(14.6587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00127419,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322468,'amu*angstrom^2'), symmetry=1, barrier=(30.7382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 8,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.132698,0.0921006,-0.000133147,9.59955e-08,-2.54016e-11,126962,34.5629], Tmin=(100,'K'), Tmax=(676.437,'K')), NASAPolynomial(coeffs=[13.2589,0.0260742,-1.24421e-05,2.37112e-09,-1.63802e-13,124921,-25.5929], Tmin=(676.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1054.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]([O])C([O])[CH][C]=C(28223)',
    structure = SMILES('[CH][C]([O])C([O])[CH][C]=C'),
    E0 = (1031.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,360,370,350,2950,3100,1380,975,1025,1650,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 8,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.118132,0.095908,-0.000137707,1.02566e-07,-3.02472e-11,124180,31.7923], Tmin=(100,'K'), Tmax=(832.269,'K')), NASAPolynomial(coeffs=[14.3122,0.0265529,-1.27067e-05,2.43685e-09,-1.69589e-13,121778,-35.1762], Tmin=(832.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1031.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(C2CsJOH) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1([O])OC1[CH][C]=C(28224)',
    structure = SMILES('[CH2]C1([O])OC1[CH][C]=C'),
    E0 = (531.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.458074,0.0863009,-9.35815e-05,5.15199e-08,-1.08119e-11,64075.1,27.9615], Tmin=(100,'K'), Tmax=(1301.62,'K')), NASAPolynomial(coeffs=[19.4146,0.0166836,-3.50441e-06,3.39244e-10,-1.27584e-14,59625.8,-70.3701], Tmin=(1301.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(C=CCJCO) + radical(CJC(O)2C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]C1OC1C(=C)[O](28112)',
    structure = SMILES('[CH2][C]C1OC1C(=C)[O]'),
    E0 = (516.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.709024,0.0799869,-8.14546e-05,4.12175e-08,-7.69142e-12,62339.2,33.9907], Tmin=(100,'K'), Tmax=(1577.55,'K')), NASAPolynomial(coeffs=[19.0666,0.012034,9.25361e-08,-4.00494e-10,3.80621e-14,58316,-63.4065], Tmin=(1577.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=C(C)OJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1([O])C([O])C1[C]=C(28202)',
    structure = SMILES('[CH2]C1([O])C([O])C1[C]=C'),
    E0 = (636.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.191178,0.07455,-6.59492e-05,2.42306e-08,-1.82325e-12,76728.9,28.8252], Tmin=(100,'K'), Tmax=(1014.42,'K')), NASAPolynomial(coeffs=[18.7239,0.0179256,-6.54808e-06,1.18078e-09,-8.28493e-14,72122.3,-65.022], Tmin=(1014.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CJC(C)2O) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])=C(O)[CH][C]=C(28225)',
    structure = SMILES('[CH2][C]=CC(O)=C([CH2])[O]'),
    E0 = (246.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.68473,0.0890838,-0.00010274,5.72762e-08,-1.19327e-11,29772.2,31.9368], Tmin=(100,'K'), Tmax=(1336.67,'K')), NASAPolynomial(coeffs=[22.1094,0.0093721,-3.82791e-07,-2.11256e-10,2.31707e-14,24705.9,-80.8031], Tmin=(1336.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])C(=O)[CH][C]=C(28226)',
    structure = SMILES('[CH2]C([O])C([O])=C[C]=C'),
    E0 = (426.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.360748,0.085591,-9.62234e-05,5.33058e-08,-1.12724e-11,51478.2,33.2266], Tmin=(100,'K'), Tmax=(1237.55,'K')), NASAPolynomial(coeffs=[20.7357,0.0131271,-3.20867e-06,4.06706e-10,-2.20898e-14,46584.1,-71.7247], Tmin=(1237.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CJCO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]([O])C(=O)C[C]=C(27787)',
    structure = SMILES('[CH2]C([O])=C([O])C[C]=C'),
    E0 = (295.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,253.166,253.413,253.52,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00263386,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356561,'amu*angstrom^2'), symmetry=1, barrier=(16.2479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356332,'amu*angstrom^2'), symmetry=1, barrier=(16.253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.237352,0.0784692,-8.5564e-05,4.50114e-08,-8.29265e-12,35737,32.4151], Tmin=(100,'K'), Tmax=(920.462,'K')), NASAPolynomial(coeffs=[16.9142,0.0176924,-5.57915e-06,8.83727e-10,-5.64768e-14,32171.5,-49.3501], Tmin=(920.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])C(=O)[CH][C]=C(28227)',
    structure = SMILES('[CH2][C]=CC([O])=C([CH2])[O]'),
    E0 = (383.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,454.926,455.733,456.852],'cm^-1')),
        HinderedRotor(inertia=(0.0725753,'amu*angstrom^2'), symmetry=1, barrier=(10.8104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0727758,'amu*angstrom^2'), symmetry=1, barrier=(10.7806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.548931,'amu*angstrom^2'), symmetry=1, barrier=(80.7729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.287324,0.0777679,-8.35762e-05,3.81496e-08,-3.70867e-12,46303.1,29.9626], Tmin=(100,'K'), Tmax=(850.559,'K')), NASAPolynomial(coeffs=[17.7438,0.0139262,-3.17759e-06,3.63001e-10,-1.78762e-14,42673.3,-55.3104], Tmin=(850.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C([O])[C]([CH2])[O](28228)',
    structure = SMILES('[CH]=C=CC([O])[C]([CH2])[O]'),
    E0 = (758.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,540,610,2055,3120,650,792.5,1650,372.869,372.943,373.044],'cm^-1')),
        HinderedRotor(inertia=(0.108046,'amu*angstrom^2'), symmetry=1, barrier=(10.6676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108118,'amu*angstrom^2'), symmetry=1, barrier=(10.6602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260225,'amu*angstrom^2'), symmetry=1, barrier=(25.7188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0568929,0.0886802,-0.000120338,8.34068e-08,-2.26062e-11,91356.5,33.2662], Tmin=(100,'K'), Tmax=(909.294,'K')), NASAPolynomial(coeffs=[15.4383,0.0210166,-8.71691e-06,1.56894e-09,-1.05566e-13,88559.3,-39.4774], Tmin=(909.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(758.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2][C]=CC=O(16857)',
    structure = SMILES('C=[C]C=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61748,'amu*angstrom^2'), symmetry=1, barrier=(37.1891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98075,0.032209,7.5369e-06,-4.54864e-08,2.38737e-11,27216,16.1108], Tmin=(100,'K'), Tmax=(899.952,'K')), NASAPolynomial(coeffs=[16.1069,0.00249496,1.93929e-06,-5.05334e-10,3.47558e-14,23334.1,-57.9917], Tmin=(899.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([O])=C[CH][C]=C(28229)',
    structure = SMILES('[CH2][C]=CC=C([CH2])[O]'),
    E0 = (459.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,432.035,432.068],'cm^-1')),
        HinderedRotor(inertia=(0.0727986,'amu*angstrom^2'), symmetry=1, barrier=(9.64662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160773,'amu*angstrom^2'), symmetry=1, barrier=(21.3233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.563729,'amu*angstrom^2'), symmetry=1, barrier=(74.661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.338117,0.0664152,-6.40179e-05,3.13594e-08,-5.85651e-12,55460.8,29.034], Tmin=(100,'K'), Tmax=(1480.54,'K')), NASAPolynomial(coeffs=[16.8531,0.0139942,-3.00326e-06,3.25919e-10,-1.52637e-14,51425.7,-54.2338], Tmin=(1480.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])[C]([O])[CH][C]=C(28230)',
    structure = SMILES('[CH2][C][CH]C([O])=C([CH2])[O]'),
    E0 = (782.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 8,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.161731,0.080082,-9.11514e-05,4.68948e-08,-7.61076e-12,94248.7,36.5951], Tmin=(100,'K'), Tmax=(887.933,'K')), NASAPolynomial(coeffs=[18.5303,0.012712,-3.31957e-06,4.54178e-10,-2.67019e-14,90380.5,-53.2529], Tmin=(887.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(CCJCO) + radical(C=C(O)CJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C([O])[C]([CH2])[O](28231)',
    structure = SMILES('[CH][C]=CC([O])[C]([CH2])[O]'),
    E0 = (1035.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 8,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.100769,0.09097,-0.000121916,8.99935e-08,-2.6889e-11,124721,35.0163], Tmin=(100,'K'), Tmax=(815.746,'K')), NASAPolynomial(coeffs=[11.9209,0.033013,-1.53497e-05,2.90708e-09,-2.01135e-13,122792,-19.6025], Tmin=(815.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1035.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]([O])[C](O)[CH][C]=C(28232)',
    structure = SMILES('[CH2][C][CH]C(O)=C([CH2])[O]'),
    E0 = (644.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.714109,0.0901604,-0.00010544,5.87686e-08,-1.22587e-11,77713.7,38.2303], Tmin=(100,'K'), Tmax=(1314.65,'K')), NASAPolynomial(coeffs=[23.1831,0.00773324,-3.04806e-07,-1.67908e-10,1.80623e-14,72270,-80.4041], Tmin=(1314.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO) + radical(C=C(O)CJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])[C]([O])[CH][C]=C(28233)',
    structure = SMILES('[CH2][C]C=C([O])C([CH2])[O]'),
    E0 = (770.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0452991,0.0827953,-9.15891e-05,5.03064e-08,-1.07198e-11,92839.7,36.9684], Tmin=(100,'K'), Tmax=(1155.43,'K')), NASAPolynomial(coeffs=[18.9255,0.0171194,-6.32649e-06,1.11062e-09,-7.52367e-14,88455.9,-57.2952], Tmin=(1155.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(770.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CJCO) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]([O])[C]([O])C[C]=C(27789)',
    structure = SMILES('[CH2][C]([O])[C]([O])C[C]=C'),
    E0 = (854.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,363.333,366.667,370,300,400,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.320106,0.10495,-0.000175523,1.52458e-07,-5.09872e-11,102905,34.3612], Tmin=(100,'K'), Tmax=(864.217,'K')), NASAPolynomial(coeffs=[11.242,0.0320272,-1.52685e-05,2.85069e-09,-1.9198e-13,101631,-15.5391], Tmin=(864.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(854.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C([O])[C]([CH2])[O](28234)',
    structure = SMILES('[CH]C=CC([O])[C]([CH2])[O]'),
    E0 = (798.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242478,0.0850858,-9.25802e-05,5.46156e-08,-1.31581e-11,96112.7,34.0341], Tmin=(100,'K'), Tmax=(998.684,'K')), NASAPolynomial(coeffs=[13.1343,0.0334511,-1.50268e-05,2.84558e-09,-1.98654e-13,93537.7,-28.1445], Tmin=(998.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(798.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C](O)[C]([O])[CH][C]=C(28235)',
    structure = SMILES('[CH2][C][CH]C([O])=C([CH2])O'),
    E0 = (644.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.714109,0.0901604,-0.00010544,5.87686e-08,-1.22587e-11,77713.7,38.2303], Tmin=(100,'K'), Tmax=(1314.65,'K')), NASAPolynomial(coeffs=[23.1831,0.00773324,-3.04806e-07,-1.67908e-10,1.80623e-14,72270,-80.4041], Tmin=(1314.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO) + radical(C=C(O)CJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=[C][CH][C]([O])[C](C)[O](28236)',
    structure = SMILES('[CH2][C][CH]C([O])=C(C)[O]'),
    E0 = (623.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.116069,0.0805552,-9.17548e-05,5.30448e-08,-1.1921e-11,75137.2,36.0658], Tmin=(100,'K'), Tmax=(1097.31,'K')), NASAPolynomial(coeffs=[17.304,0.0179,-6.10564e-06,1.00856e-09,-6.54914e-14,71365.1,-48.4519], Tmin=(1097.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(CCJCO) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]([O])[C]([O])[CH]C=C(28237)',
    structure = SMILES('[CH2][CH][CH]C([O])=C([CH2])[O]'),
    E0 = (528.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.480177,0.076418,-8.81406e-05,5.18824e-08,-1.12942e-11,63712.8,38.2839], Tmin=(100,'K'), Tmax=(857.047,'K')), NASAPolynomial(coeffs=[14.2475,0.0211636,-7.18681e-06,1.16425e-09,-7.39876e-14,61022.4,-27.94], Tmin=(857.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(C=C(O)CJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=[C]CC([O])[C]([CH2])[O](25979)',
    structure = SMILES('[CH]=[C]CC([O])[C]([CH2])[O]'),
    E0 = (924.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,331.066,331.23,331.456,331.515],'cm^-1')),
        HinderedRotor(inertia=(0.14155,'amu*angstrom^2'), symmetry=1, barrier=(11.0309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141916,'amu*angstrom^2'), symmetry=1, barrier=(11.0338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00153631,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00153449,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.2174,0.0994515,-0.000154013,1.25999e-07,-4.0536e-11,111380,34.9728], Tmin=(100,'K'), Tmax=(835.655,'K')), NASAPolynomial(coeffs=[12.769,0.0292771,-1.36668e-05,2.5594e-09,-1.74052e-13,109489,-23.6733], Tmin=(835.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(924.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(O)[C]([CH2])[O](28238)',
    structure = SMILES('[CH][C]=CC(O)[C]([CH2])[O]'),
    E0 = (805.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0951052,0.0953891,-0.000129484,9.74541e-08,-2.96596e-11,97021.7,35.9089], Tmin=(100,'K'), Tmax=(802.016,'K')), NASAPolynomial(coeffs=[12.0867,0.0346265,-1.58286e-05,2.96869e-09,-2.03934e-13,95067.9,-20.1723], Tmin=(802.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(805.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C([O])C([CH2])[O](28239)',
    structure = SMILES('[CH][C]=CC([O])C([CH2])[O]'),
    E0 = (859.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32846,0.0814289,-8.22458e-05,4.44557e-08,-9.81203e-12,103473,34.7701], Tmin=(100,'K'), Tmax=(1084.78,'K')), NASAPolynomial(coeffs=[13.655,0.0322886,-1.4296e-05,2.696e-09,-1.88e-13,100582,-30.6072], Tmin=(1084.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(859.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CJCO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C([O])[C]([CH2])O(28240)',
    structure = SMILES('[CH][C]=CC([O])[C]([CH2])O'),
    E0 = (805.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0951052,0.0953891,-0.000129484,9.74541e-08,-2.96596e-11,97021.7,35.9089], Tmin=(100,'K'), Tmax=(802.016,'K')), NASAPolynomial(coeffs=[12.0867,0.0346265,-1.58286e-05,2.96869e-09,-2.03934e-13,95067.9,-20.1723], Tmin=(802.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(805.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C([O])[C](C)[O](28241)',
    structure = SMILES('[CH][C]=CC([O])[C](C)[O]'),
    E0 = (824.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.384704,0.084148,-9.60618e-05,6.2301e-08,-1.67901e-11,99263.1,33.3689], Tmin=(100,'K'), Tmax=(890.996,'K')), NASAPolynomial(coeffs=[10.8234,0.0372865,-1.71725e-05,3.27588e-09,-2.29029e-13,97402.8,-15.7872], Tmin=(890.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4372.6,'J/mol'), sigma=(7.32671,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=682.99 K, Pc=25.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.592707,0.113023,-0.000194434,1.73225e-07,-5.92604e-11,94847.8,31.8628], Tmin=(100,'K'), Tmax=(856.735,'K')), NASAPolynomial(coeffs=[10.7935,0.0360323,-1.79144e-05,3.40066e-09,-2.31166e-13,93771.3,-16.2045], Tmin=(856.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(787.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(CJC(C)2O) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])C1C=C([CH2])O1(28242)',
    structure = SMILES('[CH2][C]([O])C1[CH]C(=C)O1'),
    E0 = (447.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35696,0.0891494,-9.18991e-05,4.48969e-08,-8.06072e-12,54049.4,28.4219], Tmin=(100,'K'), Tmax=(1621.27,'K')), NASAPolynomial(coeffs=[24.0388,0.00765036,9.36845e-07,-4.45711e-10,3.64781e-14,48291.1,-98.7323], Tmin=(1621.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1=CC([O])C1([CH2])[O](28243)',
    structure = SMILES('[CH2]C1([O])C(=C)[CH]C1[O]'),
    E0 = (511.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516271,0.066245,-4.52107e-05,9.73744e-09,1.13587e-12,61682.4,27.331], Tmin=(100,'K'), Tmax=(1101.54,'K')), NASAPolynomial(coeffs=[16.6027,0.0241319,-1.0062e-05,1.89949e-09,-1.3431e-13,57149.4,-56.3215], Tmin=(1101.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(C=CC(C)2OJ) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C=CC(=O)[C]([CH2])[O](28244)',
    structure = SMILES('[CH2]C=CC([O])=C([CH2])[O]'),
    E0 = (146.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.316255,0.0808222,-8.63307e-05,4.59028e-08,-9.21025e-12,17727.3,31.645], Tmin=(100,'K'), Tmax=(1383.13,'K')), NASAPolynomial(coeffs=[19.7559,0.012909,-1.98098e-06,8.98769e-11,2.45446e-15,13118.4,-68.2909], Tmin=(1383.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C#CC([O])[C]([CH2])[O](28245)',
    structure = SMILES('[CH2]C#CC([O])[C]([CH2])[O]'),
    E0 = (743.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2100,2250,500,550,180,754.493,754.516,3457.57],'cm^-1')),
        HinderedRotor(inertia=(0.116871,'amu*angstrom^2'), symmetry=1, barrier=(18.9196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.48341,'amu*angstrom^2'), symmetry=1, barrier=(80.0905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468423,'amu*angstrom^2'), symmetry=1, barrier=(18.9196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.48343,'amu*angstrom^2'), symmetry=1, barrier=(80.0909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282122,0.0827871,-0.000105713,6.97998e-08,-1.81195e-11,89559.7,33.0651], Tmin=(100,'K'), Tmax=(946.713,'K')), NASAPolynomial(coeffs=[14.7498,0.0216579,-8.85616e-06,1.59255e-09,-1.07564e-13,86820.4,-35.9403], Tmin=(946.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Propargyl)"""),
)

species(
    label = '[CH2]C=[C]C([O])[C]([CH2])[O](28246)',
    structure = SMILES('[CH2]C=[C]C([O])[C]([CH2])[O]'),
    E0 = (816.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180,716.555,717.519],'cm^-1')),
        HinderedRotor(inertia=(0.0304364,'amu*angstrom^2'), symmetry=1, barrier=(3.31006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.031937,'amu*angstrom^2'), symmetry=1, barrier=(11.643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0319196,'amu*angstrom^2'), symmetry=1, barrier=(11.6447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05262,'amu*angstrom^2'), symmetry=1, barrier=(24.2017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.163589,0.0877849,-0.000110932,7.32764e-08,-1.92896e-11,98358.5,33.9641], Tmin=(100,'K'), Tmax=(927.109,'K')), NASAPolynomial(coeffs=[14.2935,0.0268223,-1.23004e-05,2.35288e-09,-1.64937e-13,95738.5,-33.1355], Tmin=(927.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(816.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(O)[C]([CH2])[O](28247)',
    structure = SMILES('[CH2][C]=[C]C(O)[C]([CH2])[O]'),
    E0 = (824.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.220786,0.0987515,-0.000150683,1.20554e-07,-3.80164e-11,99269.4,35.9995], Tmin=(100,'K'), Tmax=(822.108,'K')), NASAPolynomial(coeffs=[13.5611,0.0274268,-1.27577e-05,2.39186e-09,-1.63065e-13,97147.6,-26.9131], Tmin=(822.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([O])C([CH2])[O](28248)',
    structure = SMILES('[CH2][C]=[C]C([O])C([CH2])[O]'),
    E0 = (877.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,180,387.22,738.471,738.473],'cm^-1')),
        HinderedRotor(inertia=(0.786537,'amu*angstrom^2'), symmetry=1, barrier=(18.084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786529,'amu*angstrom^2'), symmetry=1, barrier=(18.0838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00329535,'amu*angstrom^2'), symmetry=1, barrier=(37.4156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0467311,'amu*angstrom^2'), symmetry=1, barrier=(18.0839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288791,0.0836387,-9.87845e-05,6.06574e-08,-1.48618e-11,105718,34.5616], Tmin=(100,'K'), Tmax=(993.013,'K')), NASAPolynomial(coeffs=[14.5841,0.0260541,-1.17983e-05,2.25754e-09,-1.58795e-13,102879,-34.3045], Tmin=(993.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(877.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])C([O])[C]=[C]C(28249)',
    structure = SMILES('[CH2][C]([O])C([O])[C]=[C]C'),
    E0 = (903.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1670,1700,300,440,364.665,364.72,364.741,4000],'cm^-1')),
        HinderedRotor(inertia=(0.107017,'amu*angstrom^2'), symmetry=1, barrier=(10.0982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106966,'amu*angstrom^2'), symmetry=1, barrier=(10.0981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26293,'amu*angstrom^2'), symmetry=1, barrier=(24.829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106962,'amu*angstrom^2'), symmetry=1, barrier=(10.0985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0759778,0.0975099,-0.000154482,1.30838e-07,-4.35116e-11,108747,35.0051], Tmin=(100,'K'), Tmax=(835.592,'K')), NASAPolynomial(coeffs=[11.2348,0.0316301,-1.51532e-05,2.86888e-09,-1.96172e-13,107267,-15.0797], Tmin=(835.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(903.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([O])[C]([CH2])O(28250)',
    structure = SMILES('[CH2][C]=[C]C([O])[C]([CH2])O'),
    E0 = (824.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.220786,0.0987515,-0.000150683,1.20554e-07,-3.80164e-11,99269.4,35.9995], Tmin=(100,'K'), Tmax=(822.108,'K')), NASAPolynomial(coeffs=[13.5611,0.0274268,-1.27577e-05,2.39186e-09,-1.63065e-13,97147.6,-26.9131], Tmin=(822.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C([O])[C](C)[O](28251)',
    structure = SMILES('[CH2][C]=[C]C([O])[C](C)[O]'),
    E0 = (842.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1670,1700,300,440,184.67,504.902,507.167,2481.82],'cm^-1')),
        HinderedRotor(inertia=(0.0693097,'amu*angstrom^2'), symmetry=1, barrier=(12.6387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545726,'amu*angstrom^2'), symmetry=1, barrier=(12.5673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153829,'amu*angstrom^2'), symmetry=1, barrier=(28.3431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.070442,'amu*angstrom^2'), symmetry=1, barrier=(12.6145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25638,0.0874999,-0.00011702,8.47961e-08,-2.47578e-11,101511,33.4716], Tmin=(100,'K'), Tmax=(835.163,'K')), NASAPolynomial(coeffs=[12.2083,0.0302563,-1.42071e-05,2.72536e-09,-1.90431e-13,99514.6,-22.0365], Tmin=(835.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(842.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]([O])[C]([O])C=[C]C(28252)',
    structure = SMILES('[CH2]C([O])=C([O])[CH][C]C'),
    E0 = (577.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204615,0.0833155,-9.49832e-05,5.37784e-08,-1.16263e-11,69580.2,35.424], Tmin=(100,'K'), Tmax=(1217.36,'K')), NASAPolynomial(coeffs=[19.6536,0.0134383,-3.18071e-06,3.82021e-10,-1.94725e-14,65088.2,-62.8781], Tmin=(1217.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(CCJCO) + radical(C=C(O)CJ) + radical(CCJ2_triplet)"""),
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
    E0 = (794.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (794.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1319.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1253.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1241.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1266.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1243.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (797.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (797.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (797.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (817.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (817.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (817.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (794.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (985.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (896.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (794.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (794.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (994.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1247.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (908.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (952.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1012.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (989.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (905.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (952.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (919.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1070.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (923.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (892.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (923.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (956.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (958.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (802.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (802.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (858.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (970.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (933.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (1020.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (947.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (1019.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (1064.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (972.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (1011.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (925.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['C3H3(5450)', 'C=C([O])[CH][O](2850)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]([O])C([O])C=C=C(25984)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH2][C]([O])[CH][CH][C]=C(28221)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][C][O](2821)', 'C=[C][CH][CH][O](17610)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C][CH2](16918)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2][C]([O])C([O])[C][C]=C(28222)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH][C]([O])C([O])[CH][C]=C(28223)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C1([O])OC1[CH][C]=C(28224)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]C1OC1C(=C)[O](28112)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C1([O])C([O])C1[C]=C(28202)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C([O])=C(O)[CH][C]=C(28225)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C([O])C(=O)[CH][C]=C(28226)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]([O])C(=O)C[C]=C(27787)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH2][C]([O])C(=O)[CH][C]=C(28227)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(199.039,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 196.1 to 199.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C#C[CH]C([O])[C]([CH2])[O](28228)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C][O](2821)', '[CH2][C]=CC=O(16857)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.75465,'m^3/(mol*s)'), n=1.70829, Ea=(22.0129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [CO_O;YJ] for rate rule [CO_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C][CH2](16918)', 'C=C([O])[CH][O](2850)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.75465,'m^3/(mol*s)'), n=1.70829, Ea=(159.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [CO_O;YJ] for rate rule [CO_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 157.0 to 159.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(4)', '[CH2]C([O])=C[CH][C]=C(28229)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(91.7386,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 91.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2][C]([O])[C]([O])[CH][C]=C(28230)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH]=[C][CH]C([O])[C]([CH2])[O](28231)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]([O])[C](O)[CH][C]=C(28232)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C([O])[C]([O])[CH][C]=C(28233)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][C]([O])[C]([O])C[C]=C(27789)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH]=C[CH]C([O])[C]([CH2])[O](28234)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C](O)[C]([O])[CH][C]=C(28235)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.23689e+09,'s^-1'), n=1.09705, Ea=(110.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;XH_out] for rate rule [R3HJ;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['C=[C][CH][C]([O])[C](C)[O](28236)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]([O])[C]([O])[CH]C=C(28237)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.30234e+06,'s^-1'), n=1.68744, Ea=(125.264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC([O])[C]([CH2])[O](25979)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH]=[C][CH]C(O)[C]([CH2])[O](28238)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C][CH]C([O])C([CH2])[O](28239)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH]=[C][CH]C([O])[C]([CH2])O(28240)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6Hall;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=[C][CH]C([O])[C](C)[O](28241)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C([O])([CH][O])[CH][C]=C(25056)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]([O])C1C=C([CH2])O1(28242)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C1=CC([O])C1([CH2])[O](28243)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2]C=CC(=O)[C]([CH2])[O](28244)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH2]C#CC([O])[C]([CH2])[O](28245)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C3H3(5450)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C=[C]C([O])[C]([CH2])[O](28246)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][C]=[C]C(O)[C]([CH2])[O](28247)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(117344,'s^-1'), n=2.01217, Ea=(123.77,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;O_H_out] + [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][C]=[C]C([O])C([CH2])[O](28248)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][C]([O])C([O])[C]=[C]C(28249)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][C]=[C]C([O])[C]([CH2])O(28250)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1194.19,'s^-1'), n=2.56519, Ea=(148.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;O_H_out] + [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][C]=[C]C([O])[C](C)[O](28251)'],
    products = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(115.297,'s^-1'), n=2.99825, Ea=(169.04,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C]=CC([O])[C]([CH2])[O](25999)'],
    products = ['[CH2][C]([O])[C]([O])C=[C]C(28252)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #6706',
    isomers = [
        '[CH2][C]=CC([O])[C]([CH2])[O](25999)',
    ],
    reactants = [
        ('C3H3(5450)', 'C=C([O])[CH][O](2850)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #6706',
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


species(
    label = '[CH2]C([O])([CH][O])C(=C)[O](24641)',
    structure = SMILES('[CH2]C([O])([CH][O])C(=C)[O]'),
    E0 = (373.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,180,180,180,180,995.756,1242.29,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16674,'amu*angstrom^2'), symmetry=1, barrier=(3.83369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16674,'amu*angstrom^2'), symmetry=1, barrier=(3.83369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16674,'amu*angstrom^2'), symmetry=1, barrier=(3.83369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.165248,0.0924419,-0.000149674,1.2801e-07,-4.23998e-11,45090.4,34.3988], Tmin=(100,'K'), Tmax=(865.508,'K')), NASAPolynomial(coeffs=[10.4293,0.0294848,-1.36643e-05,2.52687e-09,-1.6956e-13,43895.1,-10.2782], Tmin=(865.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'ketene(1375)',
    structure = SMILES('C=C=O'),
    E0 = (-60.9858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48302,0.0083299,5.93417e-06,-1.33939e-08,5.73133e-12,-7313.53,6.51898], Tmin=(100,'K'), Tmax=(954.863,'K')), NASAPolynomial(coeffs=[5.88246,0.00580173,-1.91266e-06,3.35917e-10,-2.37398e-14,-8114.73,-6.74202], Tmin=(954.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C([O])(C=O)C(=C)[O](24653)',
    structure = SMILES('[CH2]C([O])(C=O)C(=C)[O]'),
    E0 = (76.5109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180,180,180,1600,1866.67,2653.02,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16525,'amu*angstrom^2'), symmetry=1, barrier=(3.79942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16525,'amu*angstrom^2'), symmetry=1, barrier=(3.79942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16525,'amu*angstrom^2'), symmetry=1, barrier=(3.79942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0488963,0.0984056,-0.000166867,1.44105e-07,-4.74336e-11,9339.06,30.6209], Tmin=(100,'K'), Tmax=(885.553,'K')), NASAPolynomial(coeffs=[11.2026,0.0274531,-1.25865e-05,2.28916e-09,-1.50907e-13,8135.6,-17.837], Tmin=(885.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.5109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=C(C)OJ) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C=C([O])C[C]([O])[CH][O](25554)',
    structure = SMILES('C=C([O])C[C]([O])[CH][O]'),
    E0 = (363.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,277.021,277.062,277.165,277.207,2128.7,2128.83],'cm^-1')),
        HinderedRotor(inertia=(0.00218996,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224058,'amu*angstrom^2'), symmetry=1, barrier=(12.1842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224049,'amu*angstrom^2'), symmetry=1, barrier=(12.1759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4794.68,'J/mol'), sigma=(7.69125,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=748.92 K, Pc=23.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0281748,0.0977631,-0.000167462,1.48243e-07,-5.00997e-11,43868.6,32.8408], Tmin=(100,'K'), Tmax=(873.619,'K')), NASAPolynomial(coeffs=[9.82632,0.0308335,-1.46549e-05,2.72025e-09,-1.82037e-13,42998.8,-8.28617], Tmin=(873.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]([O])C([O])C(=C)[O](25526)',
    structure = SMILES('[CH2][C]([O])C([O])C(=C)[O]'),
    E0 = (386.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,339.719,339.851,339.932,340.026,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00146077,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237682,'amu*angstrom^2'), symmetry=1, barrier=(19.4776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237906,'amu*angstrom^2'), symmetry=1, barrier=(19.4784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0334563,0.0870414,-0.000112386,7.21597e-08,-1.79869e-11,46641,33.4213], Tmin=(100,'K'), Tmax=(989.662,'K')), NASAPolynomial(coeffs=[17.2021,0.0176495,-7.21102e-06,1.31052e-09,-8.95783e-14,43242.8,-49.2291], Tmin=(989.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
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
    label = '[CH2][C]([CH][O])C(=C)[O](27023)',
    structure = SMILES('[CH2]C([O])=C([CH2])[CH][O]'),
    E0 = (354.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,455.426,455.508,455.54],'cm^-1')),
        HinderedRotor(inertia=(0.000812626,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00081245,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00081327,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83384,0.0612944,-5.8312e-05,2.82975e-08,-5.38349e-12,42765.9,28.8379], Tmin=(100,'K'), Tmax=(1285.98,'K')), NASAPolynomial(coeffs=[15.3964,0.0159991,-5.47932e-06,9.08948e-10,-5.9146e-14,39020.4,-45.0809], Tmin=(1285.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
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
    label = '[CH2][C]([O])C([O])=C[O](12602)',
    structure = SMILES('[CH2]C([O])=C([O])[CH][O]'),
    E0 = (165.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,727.222,727.227,727.227,727.231],'cm^-1')),
        HinderedRotor(inertia=(0.00831993,'amu*angstrom^2'), symmetry=1, barrier=(3.12237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135802,'amu*angstrom^2'), symmetry=1, barrier=(3.12236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05703,0.058478,-6.66932e-05,3.76854e-08,-8.14783e-12,20001,28.552], Tmin=(100,'K'), Tmax=(1201.8,'K')), NASAPolynomial(coeffs=[15.0395,0.00946163,-2.42165e-06,3.16846e-10,-1.74922e-14,16819.1,-40.7311], Tmin=(1201.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=C(O)CJ) + radical(C=CCJO)"""),
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
    label = '[CH2][C]([O])C(=C)[O](2860)',
    structure = SMILES('[CH2]C([O])=C([CH2])[O]'),
    E0 = (134.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,439.424,3272.1],'cm^-1')),
        HinderedRotor(inertia=(0.110045,'amu*angstrom^2'), symmetry=1, barrier=(15.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111135,'amu*angstrom^2'), symmetry=1, barrier=(15.1925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902014,0.0595492,-7.19852e-05,4.19793e-08,-9.05855e-12,16286.5,23.8844], Tmin=(100,'K'), Tmax=(1319.15,'K')), NASAPolynomial(coeffs=[15.5077,0.00546432,6.54129e-07,-3.60584e-10,3.24517e-14,13285.4,-47.3945], Tmin=(1319.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])([C]=C)[CH][O](27024)',
    structure = SMILES('[CH2]C([O])([C]=C)[CH][O]'),
    E0 = (688.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180,180,180,232.238,945.888,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.147625,'amu*angstrom^2'), symmetry=1, barrier=(3.39418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147625,'amu*angstrom^2'), symmetry=1, barrier=(3.39418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147625,'amu*angstrom^2'), symmetry=1, barrier=(3.39418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778114,0.0821809,-0.000144161,1.354e-07,-4.82839e-11,82902.3,31.1128], Tmin=(100,'K'), Tmax=(864.004,'K')), NASAPolynomial(coeffs=[5.71281,0.0340421,-1.66753e-05,3.14949e-09,-2.13447e-13,82993.7,13.4904], Tmin=(864.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(Cds_S)"""),
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
    label = '[CH]C([O])([CH][O])C(=C)[O](27025)',
    structure = SMILES('[CH]C([O])([CH][O])C(=C)[O]'),
    E0 = (608.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23062,0.089587,-0.000143802,1.19672e-07,-3.86246e-11,73325.9,34.2932], Tmin=(100,'K'), Tmax=(864.662,'K')), NASAPolynomial(coeffs=[11.6929,0.0249402,-1.14938e-05,2.11705e-09,-1.41703e-13,71778.2,-16.8269], Tmin=(864.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(608.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])([C][O])C(=C)[O](27026)',
    structure = SMILES('[CH2]C([O])([C][O])C(=C)[O]'),
    E0 = (654.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211684,0.0897608,-0.000143979,1.19385e-07,-3.8605e-11,78849.6,33.4154], Tmin=(100,'K'), Tmax=(843.562,'K')), NASAPolynomial(coeffs=[12.1522,0.0242111,-1.15401e-05,2.16908e-09,-1.47336e-13,77152.8,-20.2758], Tmin=(843.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CH2_triplet)"""),
)

species(
    label = 'C=C([O])C1([CH][O])CO1(27027)',
    structure = SMILES('C=C([O])C1([CH][O])CO1'),
    E0 = (120.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103555,0.0813845,-0.000103563,6.75762e-08,-1.65435e-11,14677.5,29.5582], Tmin=(100,'K'), Tmax=(1168.67,'K')), NASAPolynomial(coeffs=[15.2622,0.0160429,-2.42219e-06,2.67101e-11,1.45105e-14,12053.4,-42.0041], Tmin=(1168.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1(C(=C)[O])OC1[O](24007)',
    structure = SMILES('[CH2]C1(C(=C)[O])OC1[O]'),
    E0 = (110.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4514.26,'J/mol'), sigma=(7.30298,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=705.12 K, Pc=26.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0103366,0.0815754,-0.000101061,6.32874e-08,-1.48345e-11,13457.8,29.9318], Tmin=(100,'K'), Tmax=(1214.55,'K')), NASAPolynomial(coeffs=[16.5956,0.0139394,-1.54037e-06,-1.15912e-10,2.27916e-14,10378.9,-49.4791], Tmin=(1214.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C=C(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C([O])C1([O])CC1[O](27028)',
    structure = SMILES('C=C([O])C1([O])CC1[O]'),
    E0 = (134.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884677,0.0542852,-1.35698e-05,-3.19707e-08,1.91804e-11,16324.3,28.8535], Tmin=(100,'K'), Tmax=(942.18,'K')), NASAPolynomial(coeffs=[19.1644,0.0121545,-2.97451e-06,4.95654e-10,-3.81851e-14,11305.2,-66.6031], Tmin=(942.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CC(C)2OJ) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1([CH][O])OOC1=C(26914)',
    structure = SMILES('[CH2]C1([CH][O])OOC1=C'),
    E0 = (459.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797717,0.0647993,-5.62336e-05,2.44056e-08,-4.22926e-12,55429.1,30.3492], Tmin=(100,'K'), Tmax=(1377.63,'K')), NASAPolynomial(coeffs=[15.5853,0.0218635,-9.48454e-06,1.78292e-09,-1.23936e-13,51354.7,-45.7299], Tmin=(1377.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCOJ) + radical(CJCOOH) + radical(CCsJOH)"""),
)

species(
    label = 'C=C1OCC1([O])[CH][O](26803)',
    structure = SMILES('C=C1OCC1([O])[CH][O]'),
    E0 = (198.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.170533,0.0717229,-7.26171e-05,3.73105e-08,-7.32143e-12,24012.2,26.7824], Tmin=(100,'K'), Tmax=(1390.56,'K')), NASAPolynomial(coeffs=[17.6035,0.0144718,-3.19635e-06,3.5432e-10,-1.67445e-14,19850.8,-60.5993], Tmin=(1390.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])C(=C)OC1[O](27007)',
    structure = SMILES('[CH2]C1([O])C(=C)OC1[O]'),
    E0 = (191.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849658,0.0584483,-2.87628e-05,-1.40745e-08,1.27036e-11,23114.2,25.4347], Tmin=(100,'K'), Tmax=(926.247,'K')), NASAPolynomial(coeffs=[17.3295,0.0151999,-3.93878e-06,6.01324e-10,-4.106e-14,18863.6,-59.2733], Tmin=(926.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C([O])([CH][O])[C]1CO1(27029)',
    structure = SMILES('[CH2]C([O])([CH][O])[C]1CO1'),
    E0 = (525.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.223241,0.105816,-0.000188079,1.65641e-07,-5.37737e-11,63375.3,30.6584], Tmin=(100,'K'), Tmax=(946.942,'K')), NASAPolynomial(coeffs=[9.4357,0.029644,-1.13892e-05,1.80201e-09,-1.04713e-13,63131.9,-7.04011], Tmin=(946.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C2CsJO) + radical(CJC(C)2O) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([CH][O])OC[C]1[O](27030)',
    structure = SMILES('[CH2]C1([CH][O])OC[C]1[O]'),
    E0 = (516.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.125604,0.0983531,-0.000159442,1.30519e-07,-3.98785e-11,62285.9,28.6454], Tmin=(100,'K'), Tmax=(974.834,'K')), NASAPolynomial(coeffs=[11.8329,0.0246659,-8.17765e-06,1.16718e-09,-6.1983e-14,61124.1,-22.7432], Tmin=(974.834,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1([O])CC[C]1[O](27031)',
    structure = SMILES('[O][CH]C1([O])CC[C]1[O]'),
    E0 = (509.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807354,0.0751333,-9.112e-05,6.19798e-08,-1.7386e-11,61359.1,27.7823], Tmin=(100,'K'), Tmax=(859.696,'K')), NASAPolynomial(coeffs=[10.1921,0.0314675,-1.49309e-05,2.89691e-09,-2.04462e-13,59745.5,-16.0747], Tmin=(859.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])[C]([O])CC1[O](27032)',
    structure = SMILES('[CH2]C1([O])[C]([O])CC1[O]'),
    E0 = (532.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255552,0.0773957,-8.18578e-05,4.31318e-08,-8.89725e-12,64241.7,27.9466], Tmin=(100,'K'), Tmax=(1184.96,'K')), NASAPolynomial(coeffs=[17.599,0.0188492,-7.74412e-06,1.43401e-09,-9.97678e-14,60131.5,-58.6685], Tmin=(1184.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])[CH]OC[C]1[O](27033)',
    structure = SMILES('[CH2]C1([O])[CH]OC[C]1[O]'),
    E0 = (443.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280742,0.0788365,-8.01866e-05,2.85627e-08,2.72156e-12,53431,24.19], Tmin=(100,'K'), Tmax=(794.442,'K')), NASAPolynomial(coeffs=[17.8544,0.0138311,-1.77705e-06,-3.69774e-11,1.55272e-14,49897.9,-61.2116], Tmin=(794.442,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOCs) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])OC1([CH2])[CH][O](26930)',
    structure = SMILES('[CH2]C1([O])OC1([CH2])[CH][O]'),
    E0 = (511.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.629269,0.108306,-0.000175736,1.39973e-07,-4.16478e-11,61652.8,28.8218], Tmin=(100,'K'), Tmax=(980.295,'K')), NASAPolynomial(coeffs=[15.1449,0.0208383,-6.54607e-06,8.72015e-10,-4.23287e-14,59670.2,-41.3037], Tmin=(980.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOH) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1([O])CC1([O])[CH][O](26825)',
    structure = SMILES('[CH2]C1([O])CC1([O])[CH][O]'),
    E0 = (529.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0223684,0.0925677,-0.000134524,1.02266e-07,-3.07742e-11,63838,27.6263], Tmin=(100,'K'), Tmax=(816.349,'K')), NASAPolynomial(coeffs=[13.6197,0.0259425,-1.2104e-05,2.29234e-09,-1.5811e-13,61618,-35.214], Tmin=(816.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])C([O])C1([CH2])[O](27034)',
    structure = SMILES('[CH2]C1([O])C([O])C1([CH2])[O]'),
    E0 = (553.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.32993,0.0923603,-0.000116113,7.0889e-08,-1.6665e-11,66712.2,27.7764], Tmin=(100,'K'), Tmax=(1052.01,'K')), NASAPolynomial(coeffs=[19.9912,0.0150937,-5.94203e-06,1.07232e-09,-7.36162e-14,62436.6,-71.2916], Tmin=(1052.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CC(C)2OJ) + radical(CC(C)OJ) + radical(CJC(C)2O) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])[CH]OC1([CH2])[O](26868)',
    structure = SMILES('[CH2]C1([O])[CH]OC1([CH2])[O]'),
    E0 = (504.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.32937,0.0906359,-9.89388e-05,3.67102e-08,2.75603e-12,60867.5,26.8019], Tmin=(100,'K'), Tmax=(803.846,'K')), NASAPolynomial(coeffs=[22.4987,0.00683759,1.83148e-06,-7.51703e-10,6.59299e-14,56234.8,-84.3336], Tmin=(803.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)2OJ) + radical(CC(C)(O)OJ) + radical(CCsJOCs) + radical(CJC(C)2O) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C([O])C([O])=C[O](11158)',
    structure = SMILES('C=C([O])C([O])=C[O]'),
    E0 = (-125.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04411,'amu*angstrom^2'), symmetry=1, barrier=(24.0062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4714.77,'J/mol'), sigma=(7.20548,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=736.44 K, Pc=28.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0729095,0.0724912,-8.7986e-05,4.83801e-08,-9.6461e-12,-14943,25.3649], Tmin=(100,'K'), Tmax=(1456.33,'K')), NASAPolynomial(coeffs=[21.1206,-0.00117263,3.8032e-06,-9.24506e-10,6.85246e-14,-19477.2,-79.2224], Tmin=(1456.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C([O])C(=C)[O](2342)',
    structure = SMILES('C=C([O])C(=C)[O]'),
    E0 = (-58.255,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.748412,'amu*angstrom^2'), symmetry=1, barrier=(17.2075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4088.81,'J/mol'), sigma=(6.45843,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=638.66 K, Pc=34.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997617,0.0584393,-6.83658e-05,3.90433e-08,-8.38052e-12,-6891.8,19.5611], Tmin=(100,'K'), Tmax=(1280.26,'K')), NASAPolynomial(coeffs=[15.4008,0.00708459,-7.52445e-07,-4.1341e-11,8.63455e-15,-10059.1,-51.4511], Tmin=(1280.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(=C)[CH][O](22500)',
    structure = SMILES('[CH2]C(=C[O])C(=C)[O]'),
    E0 = (64.4117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30147,'amu*angstrom^2'), symmetry=1, barrier=(29.9235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29982,'amu*angstrom^2'), symmetry=1, barrier=(29.8854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.70548,0.0763064,-8.01811e-05,3.90501e-08,-6.92726e-12,7939.82,27.8012], Tmin=(100,'K'), Tmax=(1656,'K')), NASAPolynomial(coeffs=[21.9763,0.00356412,1.97243e-06,-5.70396e-10,4.25168e-14,2889.6,-85.6326], Tmin=(1656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=O(1376)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,623.763,623.847],'cm^-1')),
        HinderedRotor(inertia=(0.000434039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44947,0.00889162,4.91274e-06,-1.14995e-08,4.62261e-12,19378.1,9.99239], Tmin=(100,'K'), Tmax=(1030.77,'K')), NASAPolynomial(coeffs=[6.10651,0.00625992,-2.43247e-06,4.78594e-10,-3.54632e-14,18422.4,-4.88573], Tmin=(1030.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJC=O) + radical(CsCJ=O)"""),
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
    label = '[CH]=C([O])C([CH2])([O])[CH][O](25801)',
    structure = SMILES('[CH]=C([O])C([CH2])([O])[CH][O]'),
    E0 = (620.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,418.951,748.418,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156207,'amu*angstrom^2'), symmetry=1, barrier=(3.59151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156207,'amu*angstrom^2'), symmetry=1, barrier=(3.59151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156207,'amu*angstrom^2'), symmetry=1, barrier=(3.59151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0741845,0.0961701,-0.000167045,1.46278e-07,-4.86602e-11,74810.9,35.1286], Tmin=(100,'K'), Tmax=(882.327,'K')), NASAPolynomial(coeffs=[10.7731,0.0264732,-1.25261e-05,2.3033e-09,-1.52624e-13,73747.9,-10.4731], Tmin=(882.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C([CH2])([O])[CH][O](27035)',
    structure = SMILES('[CH]=C(O)C([CH2])([O])[CH][O]'),
    E0 = (483.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,208.45,1546.27,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.144841,'amu*angstrom^2'), symmetry=1, barrier=(3.33019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144841,'amu*angstrom^2'), symmetry=1, barrier=(3.33019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144841,'amu*angstrom^2'), symmetry=1, barrier=(3.33019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144841,'amu*angstrom^2'), symmetry=1, barrier=(3.33019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.269279,0.099852,-0.000158379,1.27475e-07,-3.97858e-11,58252.9,34.8617], Tmin=(100,'K'), Tmax=(862.729,'K')), NASAPolynomial(coeffs=[14.4391,0.0232408,-1.05439e-05,1.92922e-09,-1.2867e-13,56028.2,-32.111], Tmin=(862.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])(O)[CH][O](27036)',
    structure = SMILES('[CH]=C([O])C([CH2])(O)[CH][O]'),
    E0 = (391.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,180,1600,1700.59,2818.51,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153563,'amu*angstrom^2'), symmetry=1, barrier=(3.53072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.32095,0.102617,-0.000170713,1.41176e-07,-4.44121e-11,47272.6,35.3332], Tmin=(100,'K'), Tmax=(904.386,'K')), NASAPolynomial(coeffs=[13.9402,0.0225387,-9.69423e-06,1.69317e-09,-1.08198e-13,45388.5,-28.191], Tmin=(904.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CC(C)(O)CJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(C)([O])[CH][O](27037)',
    structure = SMILES('[CH]=C([O])C(C)([O])[CH][O]'),
    E0 = (407.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,180,180,180,557.472,601.889,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162585,'amu*angstrom^2'), symmetry=1, barrier=(3.73814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162585,'amu*angstrom^2'), symmetry=1, barrier=(3.73814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162585,'amu*angstrom^2'), symmetry=1, barrier=(3.73814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215433,0.0899807,-0.000140689,1.17067e-07,-3.80489e-11,49137.4,33.7469], Tmin=(100,'K'), Tmax=(858.793,'K')), NASAPolynomial(coeffs=[11.0326,0.0281946,-1.28546e-05,2.3703e-09,-1.59262e-13,47699.9,-14.345], Tmin=(858.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])([O])C[O](27038)',
    structure = SMILES('[CH]=C([O])C([CH2])([O])C[O]'),
    E0 = (440.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,180,180,509.058,655.452,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160476,'amu*angstrom^2'), symmetry=1, barrier=(3.68965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160476,'amu*angstrom^2'), symmetry=1, barrier=(3.68965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160476,'amu*angstrom^2'), symmetry=1, barrier=(3.68965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335972,0.0875194,-0.000136175,1.14517e-07,-3.77689e-11,53119.4,34.6401], Tmin=(100,'K'), Tmax=(851.045,'K')), NASAPolynomial(coeffs=[10.1888,0.0296127,-1.36717e-05,2.54187e-09,-1.71866e-13,51862.3,-8.83773], Tmin=(851.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])([O])[CH]O(27039)',
    structure = SMILES('[CH]=C([O])C([CH2])([O])[CH]O'),
    E0 = (395.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,180,180,180,1600,1659.18,2852.93,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151016,'amu*angstrom^2'), symmetry=1, barrier=(3.47215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151016,'amu*angstrom^2'), symmetry=1, barrier=(3.47215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151016,'amu*angstrom^2'), symmetry=1, barrier=(3.47215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151016,'amu*angstrom^2'), symmetry=1, barrier=(3.47215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.282609,0.0988497,-0.000154496,1.2082e-07,-3.63695e-11,47682.6,35.5489], Tmin=(100,'K'), Tmax=(894.938,'K')), NASAPolynomial(coeffs=[15.4681,0.0200487,-8.33637e-06,1.45067e-09,-9.33776e-14,45199.9,-36.8112], Tmin=(894.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])([CH][O])C[C]=O(24559)',
    structure = SMILES('[CH2]C([O])([CH][O])C[C]=O'),
    E0 = (405.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,180,180,452.256,719.631,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157547,'amu*angstrom^2'), symmetry=1, barrier=(3.62232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157547,'amu*angstrom^2'), symmetry=1, barrier=(3.62232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157547,'amu*angstrom^2'), symmetry=1, barrier=(3.62232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157547,'amu*angstrom^2'), symmetry=1, barrier=(3.62232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4579.1,'J/mol'), sigma=(7.44027,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=715.24 K, Pc=25.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.28415,0.111412,-0.000214242,2.03029e-07,-7.1372e-11,48936,33.2145], Tmin=(100,'K'), Tmax=(882.865,'K')), NASAPolynomial(coeffs=[7.39975,0.0363446,-1.831e-05,3.43402e-09,-2.29113e-13,49148.1,5.986], Tmin=(882.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CJC(C)2O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([O])([C]=O)[CH][O](27040)',
    structure = SMILES('[CH2]C([O])([C]=O)[CH][O]'),
    E0 = (448.989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1855,455,950,180,180,180,180,1600,1753.38,2778.86,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157181,'amu*angstrom^2'), symmetry=1, barrier=(3.6139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157181,'amu*angstrom^2'), symmetry=1, barrier=(3.6139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157181,'amu*angstrom^2'), symmetry=1, barrier=(3.6139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267849,0.0999325,-0.000208605,2.00797e-07,-7.02322e-11,54118.5,30.0536], Tmin=(100,'K'), Tmax=(899.518,'K')), NASAPolynomial(coeffs=[6.94477,0.0266676,-1.377e-05,2.54466e-09,-1.65605e-13,54680.1,8.34721], Tmin=(899.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)(C=O)OJ) + radical(CCOJ) + radical(CJC(C)(C=O)O) + radical(CCsJOH) + radical(CC(C)(O)CJ=O)"""),
)

species(
    label = '[CH2]C1([CH][O])OCC1=O(25541)',
    structure = SMILES('[CH2]C1([CH][O])OCC1=O'),
    E0 = (198.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521158,0.0807688,-0.00010288,6.99857e-08,-1.91597e-11,23961,25.9017], Tmin=(100,'K'), Tmax=(888.906,'K')), NASAPolynomial(coeffs=[12.3719,0.0274389,-1.28833e-05,2.48602e-09,-1.74883e-13,21854.2,-29.8751], Tmin=(888.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1([O])CCC1=O(25557)',
    structure = SMILES('[O][CH]C1([O])CCC1=O'),
    E0 = (168.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57508,0.0601657,-5.55287e-05,1.6204e-08,8.8356e-12,20339.4,27.1047], Tmin=(100,'K'), Tmax=(562.997,'K')), NASAPolynomial(coeffs=[6.64051,0.0347895,-1.61946e-05,3.10965e-09,-2.17905e-13,19600.8,4.08316], Tmin=(562.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CC(C)(C=O)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])C(=O)CC1[O](26979)',
    structure = SMILES('[CH2]C1([O])C(=O)CC1[O]'),
    E0 = (191.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420408,0.0692894,-6.74989e-05,3.3049e-08,-6.30738e-12,23200.9,29.165], Tmin=(100,'K'), Tmax=(1285.19,'K')), NASAPolynomial(coeffs=[17.4588,0.0162601,-5.60668e-06,9.43982e-10,-6.22455e-14,18821.4,-57.3103], Tmin=(1285.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CC(C)(C=O)OJ) + radical(CC(C)OJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2][C]1OOC1([CH2])[CH][O](27041)',
    structure = SMILES('[CH2][C]1OOC1([CH2])[CH][O]'),
    E0 = (730.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.196257,0.0968707,-0.000171974,1.61061e-07,-5.72001e-11,87971.5,29.7512], Tmin=(100,'K'), Tmax=(863.805,'K')), NASAPolynomial(coeffs=[6.56967,0.0380244,-1.88505e-05,3.57077e-09,-2.42165e-13,87964.8,6.27083], Tmin=(863.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(730.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(CCOJ) + radical(C2CsJOO) + radical(CJCOOH) + radical(CCsJOH) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][C]1OCC1([O])[CH][O](27042)',
    structure = SMILES('[CH2][C]1OCC1([O])[CH][O]'),
    E0 = (519.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0197766,0.099913,-0.000175194,1.54489e-07,-5.0217e-11,62592.7,28.9002], Tmin=(100,'K'), Tmax=(950.387,'K')), NASAPolynomial(coeffs=[8.4446,0.030301,-1.142e-05,1.794e-09,-1.03958e-13,62533.8,-3.20122], Tmin=(950.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C2CsJOCs) + radical(CCsJOH) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1OC([O])C1([CH2])[O](27043)',
    structure = SMILES('[CH2][C]1OC([O])C1([CH2])[O]'),
    E0 = (510.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0469253,0.0967224,-0.000161631,1.38248e-07,-4.42067e-11,61549.5,28.8554], Tmin=(100,'K'), Tmax=(944.535,'K')), NASAPolynomial(coeffs=[9.75774,0.0287161,-1.09408e-05,1.75592e-09,-1.04597e-13,60914.2,-11.0918], Tmin=(944.535,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(C2CsJOCs) + radical(CJC(C)2O) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1OO[CH]C1([CH2])[O](27044)',
    structure = SMILES('[CH2][C]1OO[CH]C1([CH2])[O]'),
    E0 = (652.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440005,0.0772324,-8.16129e-05,3.68447e-08,-2.46951e-12,78580.4,27.4408], Tmin=(100,'K'), Tmax=(787.617,'K')), NASAPolynomial(coeffs=[15.3022,0.0194339,-5.2094e-06,6.75632e-10,-3.57315e-14,75690.9,-44.1937], Tmin=(787.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CC(C)2OJ) + radical(C2CsJOOC) + radical(CCsJOOC) + radical(CJC(C)2O) + radical(CJCOOH)"""),
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
    E0 = (373.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (373.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (531.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (549.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (873.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (603.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (620.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1207.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (820.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (866.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (379.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (376.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (379.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (459.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (382.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (382.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (605.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (517.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (509.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (532.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (448.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (513.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (530.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (553.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (505.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (373.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (394.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (373.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (373.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (556.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (732.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (832.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (675.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (436.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (451.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (484.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (428.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (651.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (886.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (382.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (382.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (382.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (730.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (520.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (615.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (652.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['ketene(1375)', 'C=C([O])[CH][O](2850)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C([O])(C=O)C(=C)[O](24653)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['C=C([O])C[C]([O])[CH][O](25554)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][C]([O])C([O])C(=C)[O](25526)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2][C]([CH][O])C(=C)[O](27023)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][C]([O])C([O])=C[O](12602)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][O](1548)', '[CH2][C]([O])C(=C)[O](2860)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(4)', '[CH2]C([O])([C]=C)[CH][O](27024)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([O])([CH][O])C(=C)[O](27025)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2]C([O])([C][O])C(=C)[O](27026)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['C=C([O])C1([CH][O])CO1(27027)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1(C(=C)[O])OC1[O](24007)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['C=C([O])C1([O])CC1[O](27028)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.01367e+13,'s^-1'), n=0.0238333, Ea=(5.50893,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_2H;Ypri_rad_out] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([CH][O])OOC1=C(26914)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(86.0479,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination
Ea raised from 82.6 to 86.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['C=C1OCC1([O])[CH][O](26803)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])C(=C)OC1[O](27007)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C([O])([CH][O])[C]1CO1(27029)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([CH][O])OC[C]1[O](27030)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(143.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[O][CH]C1([O])CC[C]1[O](27031)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(135.429,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 133.8 to 135.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])[C]([O])CC1[O](27032)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(159.16,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 157.6 to 159.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])[CH]OC[C]1[O](27033)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])OC1([CH2])[CH][O](26930)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(139.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])CC1([O])[CH][O](26825)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(156.303,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])C([O])C1([CH2])[O](27034)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(180.142,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])[CH]OC1([CH2])[O](26868)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(131.687,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(T)(20)', 'C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(116.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;Y_1centerbirad] for rate rule [CO_O;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.8 to 117.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH][O](1548)', 'C=C([O])C(=C)[O](2342)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(15.5093,'m^3/(mol*s)'), n=1.70829, Ea=(22.0129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [CO_O;YJ] for rate rule [CO_O;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)', 'C=C([O])C(=C)[CH][O](22500)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(66.4028,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 66.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=O(1376)', 'C=C([O])[CH][O](2850)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(191.932,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 189.3 to 191.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['ketene(1375)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=O(1376)', '[CH2][C]([O])[CH][O](26240)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH]=C([O])C([CH2])([O])[CH][O](25801)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C(O)C([CH2])([O])[CH][O](27035)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C([O])C([CH2])(O)[CH][O](27036)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C([O])C(C)([O])[CH][O](27037)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([O])C([CH2])([O])C[O](27038)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C([O])C([CH2])([O])[CH]O(27039)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C([O])([CH][O])C[C]=O(24559)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(T)(20)', '[CH2]C([O])([C]=O)[CH][O](27040)'],
    products = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([CH][O])OCC1=O(25541)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[O][CH]C1([O])CCC1=O(25557)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2]C1([O])C(=O)CC1[O](26979)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2][C]1OOC1([CH2])[CH][O](27041)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(357.094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2][C]1OCC1([O])[CH][O](27042)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(146.452,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2][C]1OC([O])C1([CH2])[O](27043)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([O])([CH][O])C(=C)[O](24641)'],
    products = ['[CH2][C]1OO[CH]C1([CH2])[O](27044)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6.44463e+09,'s^-1'), n=0.470283, Ea=(278.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 277.2 to 278.5 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #6165',
    isomers = [
        '[CH2]C([O])([CH][O])C(=C)[O](24641)',
    ],
    reactants = [
        ('ketene(1375)', 'C=C([O])[CH][O](2850)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #6165',
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


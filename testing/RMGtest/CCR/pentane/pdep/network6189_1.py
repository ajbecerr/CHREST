species(
    label = 'C#COOC([CH2])[CH][O](25545)',
    structure = SMILES('C#COOC([CH2])[CH][O]'),
    E0 = (575.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180,2243.53],'cm^-1')),
        HinderedRotor(inertia=(0.426618,'amu*angstrom^2'), symmetry=1, barrier=(9.80878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60797,'amu*angstrom^2'), symmetry=1, barrier=(36.9704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60805,'amu*angstrom^2'), symmetry=1, barrier=(36.9721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6086,'amu*angstrom^2'), symmetry=1, barrier=(36.9849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60858,'amu*angstrom^2'), symmetry=1, barrier=(36.9845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.148531,0.10212,-0.000174214,1.56756e-07,-5.47186e-11,69350.7,32.8526], Tmin=(100,'K'), Tmax=(833.728,'K')), NASAPolynomial(coeffs=[9.69425,0.0351324,-1.81355e-05,3.51861e-09,-2.4317e-13,68396.4,-8.72374], Tmin=(833.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CJCOOH) + radical(CCsJOH)"""),
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
    label = 'C#COOC([CH2])C=O(28447)',
    structure = SMILES('C#COOC([CH2])C=O'),
    E0 = (249.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.598471,0.0721309,-7.00913e-05,3.36375e-08,-6.42262e-12,30081.2,31.5495], Tmin=(100,'K'), Tmax=(1257.76,'K')), NASAPolynomial(coeffs=[16.3024,0.0221884,-1.05301e-05,2.0675e-09,-1.4758e-13,26130.9,-47.8143], Tmin=(1257.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCOOH)"""),
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
    label = '[CH]=[C]OOC=C[O](11850)',
    structure = SMILES('[CH]=[C]OO[CH]C=O'),
    E0 = (478.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0772868,'amu*angstrom^2'), symmetry=1, barrier=(1.77698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59592,'amu*angstrom^2'), symmetry=1, barrier=(36.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.06212,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29353,0.0648156,-9.63724e-05,8.10379e-08,-2.78527e-11,57698.5,28.5617], Tmin=(100,'K'), Tmax=(744.305,'K')), NASAPolynomial(coeffs=[8.15805,0.0257386,-1.32151e-05,2.60878e-09,-1.84242e-13,56737.2,-2.1218], Tmin=(744.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OOC=C(2667)',
    structure = SMILES('[CH]=[C]OOC=C'),
    E0 = (549.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.719942,'amu*angstrom^2'), symmetry=1, barrier=(16.5529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.721218,'amu*angstrom^2'), symmetry=1, barrier=(16.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9185,'amu*angstrom^2'), symmetry=1, barrier=(21.1181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73695,0.0485798,-4.82298e-05,2.46636e-08,-5.05782e-12,66230.5,25.6851], Tmin=(100,'K'), Tmax=(1173.78,'K')), NASAPolynomial(coeffs=[11.1557,0.0164823,-7.2115e-06,1.36645e-09,-9.58103e-14,64019.4,-21.2642], Tmin=(1173.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = 'C#COOC([CH2])[C][O](28448)',
    structure = SMILES('C#COOC([CH2])[C][O]'),
    E0 = (856.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0995002,0.0994142,-0.000168481,1.48194e-07,-5.10269e-11,103110,31.8596], Tmin=(100,'K'), Tmax=(809.655,'K')), NASAPolynomial(coeffs=[11.3749,0.0299331,-1.60554e-05,3.17143e-09,-2.21837e-13,101671,-18.4852], Tmin=(809.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(856.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CJCOOH) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]C([CH][O])OOC#C(28449)',
    structure = SMILES('[CH]C([CH][O])OOC#C'),
    E0 = (809.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.135985,0.102463,-0.0001806,1.63887e-07,-5.72158e-11,97523.7,31.8162], Tmin=(100,'K'), Tmax=(839.945,'K')), NASAPolynomial(coeffs=[9.89189,0.0328792,-1.73545e-05,3.378e-09,-2.33198e-13,96609.2,-10.2297], Tmin=(839.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(809.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C#COOC1CC1[O](28450)',
    structure = SMILES('C#COOC1CC1[O]'),
    E0 = (335.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.408339,0.0683262,-5.368e-05,1.50813e-08,1.97149e-13,40526.8,26.7772], Tmin=(100,'K'), Tmax=(1075.81,'K')), NASAPolynomial(coeffs=[18.5551,0.018102,-7.70081e-06,1.49102e-09,-1.07799e-13,35624.3,-66.7351], Tmin=(1075.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Ct-CtOs) + group(Ct-CtH) + ring(Cyclopropane) + radical(CC(C)OJ)"""),
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
    label = '[O][CH]C1CO1(28451)',
    structure = SMILES('[O][CH]C1CO1'),
    E0 = (142.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15395,0.0412018,-4.99219e-05,3.58319e-08,-9.82116e-12,17234.2,16.4918], Tmin=(100,'K'), Tmax=(1096.88,'K')), NASAPolynomial(coeffs=[6.20023,0.0172046,-4.46729e-06,5.24037e-10,-2.31476e-14,16902.5,-0.869076], Tmin=(1096.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1OC1[O](23880)',
    structure = SMILES('[CH2]C1OC1[O]'),
    E0 = (133.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04577,0.036085,-3.15326e-05,1.60953e-08,-3.19114e-12,16148.6,17.7275], Tmin=(100,'K'), Tmax=(1460.65,'K')), NASAPolynomial(coeffs=[7.90194,0.0141418,-2.93314e-06,2.73796e-10,-9.37935e-15,15067.8,-10.5873], Tmin=(1460.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = 'C#COOC(C)=C[O](28452)',
    structure = SMILES('C#COOC(C)=C[O]'),
    E0 = (185.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294138,0.077679,-8.26609e-05,4.34393e-08,-8.95446e-12,22452.1,28.2148], Tmin=(100,'K'), Tmax=(1182.63,'K')), NASAPolynomial(coeffs=[17.5852,0.0191954,-8.48239e-06,1.62359e-09,-1.1487e-13,18362.3,-58.1049], Tmin=(1182.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = 'C#COOC(=C)C[O](28453)',
    structure = SMILES('C#COOC(=C)C[O]'),
    E0 = (325.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70937,0.0810833,-0.000126734,1.16161e-07,-4.26527e-11,39251.4,30.0655], Tmin=(100,'K'), Tmax=(797.157,'K')), NASAPolynomial(coeffs=[6.242,0.0387939,-1.98226e-05,3.8882e-09,-2.7259e-13,38830.9,7.52322], Tmin=(797.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ)"""),
)

species(
    label = '[O][CH]C1CC=[C]OO1(28454)',
    structure = SMILES('[O][CH]C1CC=[C]OO1'),
    E0 = (411.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42725,0.0586417,-5.55242e-05,3.20039e-08,-7.94221e-12,49608.8,25.7317], Tmin=(100,'K'), Tmax=(950.054,'K')), NASAPolynomial(coeffs=[7.76017,0.0319781,-1.34262e-05,2.4629e-09,-1.68672e-13,48405.5,-4.4964], Tmin=(950.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(CCOJ) + radical(CCsJOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1OO[C]=CC1[O](28455)',
    structure = SMILES('[CH2]C1OO[C]=CC1[O]'),
    E0 = (437.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21384,0.0492354,-7.98625e-06,-3.09535e-08,1.73664e-11,52677.1,27.7531], Tmin=(100,'K'), Tmax=(941.505,'K')), NASAPolynomial(coeffs=[15.7384,0.0174994,-5.17549e-06,8.67959e-10,-6.13663e-14,48613.7,-48.4985], Tmin=(941.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(CC(C)OJ) + radical(CJCOOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1[CH]OC=[C]OO1(28456)',
    structure = SMILES('[CH2]C1[CH]OC=[C]OO1'),
    E0 = (431.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404751,0.0212866,0.000177962,-3.06885e-07,1.37556e-10,52082.6,30.1539], Tmin=(100,'K'), Tmax=(896.967,'K')), NASAPolynomial(coeffs=[49.6554,-0.0444622,3.05749e-05,-6.07396e-09,4.05977e-13,37057,-236.604], Tmin=(896.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCsJOC(O)) + radical(CJCOOH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1CC([CH][O])OO1(28457)',
    structure = SMILES('[CH]=C1CC([CH][O])OO1'),
    E0 = (436.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12067,0.0559129,-3.95438e-05,1.36109e-08,-1.87335e-12,52621.8,29.3066], Tmin=(100,'K'), Tmax=(1693.61,'K')), NASAPolynomial(coeffs=[15.5941,0.0217294,-9.2682e-06,1.69336e-09,-1.14165e-13,47719.3,-48.1449], Tmin=(1693.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OOC([CH2])C1[O](28458)',
    structure = SMILES('[CH]=C1OOC([CH2])C1[O]'),
    E0 = (461.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22474,0.0429669,1.93301e-05,-6.24059e-08,2.82577e-11,55676.1,30.1764], Tmin=(100,'K'), Tmax=(971.035,'K')), NASAPolynomial(coeffs=[18.7299,0.0148324,-5.13926e-06,1.03094e-09,-8.17391e-14,50203.3,-64.4368], Tmin=(971.035,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1O[CH]C([CH2])OO1(28459)',
    structure = SMILES('[CH]=C1O[CH]C([CH2])OO1'),
    E0 = (465.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146428,0.0689598,-4.01617e-05,-6.06857e-09,8.44155e-12,56118.8,25.0666], Tmin=(100,'K'), Tmax=(1054.27,'K')), NASAPolynomial(coeffs=[22.9234,0.0141832,-7.24486e-06,1.58374e-09,-1.23457e-13,49557.8,-94.3624], Tmin=(1054.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(CCsJOC(O)) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = 'C#COOC(=C)[CH][O](24131)',
    structure = SMILES('C#COOC([CH2])=C[O]'),
    E0 = (344.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.30452,'amu*angstrom^2'), symmetry=1, barrier=(29.9936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30466,'amu*angstrom^2'), symmetry=1, barrier=(29.9967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30427,'amu*angstrom^2'), symmetry=1, barrier=(29.9877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30457,'amu*angstrom^2'), symmetry=1, barrier=(29.9946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0951279,0.0801655,-9.27855e-05,5.18071e-08,-1.11372e-11,41574.2,29.6172], Tmin=(100,'K'), Tmax=(1148.58,'K')), NASAPolynomial(coeffs=[19.3653,0.0130566,-5.14458e-06,9.3833e-10,-6.51931e-14,37147.5,-66.0198], Tmin=(1148.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C#COOC=C(2350)',
    structure = SMILES('C#COOC=C'),
    E0 = (294.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.4496,'amu*angstrom^2'), symmetry=1, barrier=(33.3292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46336,'amu*angstrom^2'), symmetry=1, barrier=(33.6454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45788,'amu*angstrom^2'), symmetry=1, barrier=(33.5196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.43,'J/mol'), sigma=(5.36906,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.96 K, Pc=47.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32803,0.0526921,-5.09704e-05,2.41725e-08,-4.47989e-12,35540.9,21.2321], Tmin=(100,'K'), Tmax=(1314.74,'K')), NASAPolynomial(coeffs=[14.5138,0.0125762,-5.20248e-06,9.65302e-10,-6.70875e-14,32073.7,-45.9899], Tmin=(1314.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH)"""),
)

species(
    label = 'C#COOC=C[O](11169)',
    structure = SMILES('C#COO[CH]C=O'),
    E0 = (223.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.77237,'amu*angstrom^2'), symmetry=1, barrier=(40.7503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77226,'amu*angstrom^2'), symmetry=1, barrier=(40.7478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77123,'amu*angstrom^2'), symmetry=1, barrier=(40.7241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7714,'amu*angstrom^2'), symmetry=1, barrier=(40.7279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3598.19,'J/mol'), sigma=(5.81028,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.03 K, Pc=41.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34612,0.0632936,-7.83958e-05,5.1921e-08,-1.40909e-11,26989.1,22.465], Tmin=(100,'K'), Tmax=(887.497,'K')), NASAPolynomial(coeffs=[9.9741,0.0244068,-1.26716e-05,2.55062e-09,-1.8377e-13,25457.6,-18.1302], Tmin=(887.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(OCJC=O)"""),
)

species(
    label = 'C#CO[O](2880)',
    structure = SMILES('C#CO[O]'),
    E0 = (341.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,207.446,207.513],'cm^-1')),
        HinderedRotor(inertia=(1.42756,'amu*angstrom^2'), symmetry=1, barrier=(43.6063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84093,0.0270252,-4.16758e-05,3.32096e-08,-1.02948e-11,41097.9,11.7896], Tmin=(100,'K'), Tmax=(872.981,'K')), NASAPolynomial(coeffs=[6.69652,0.00695047,-3.0441e-06,5.47527e-10,-3.61761e-14,40516.5,-5.76197], Tmin=(872.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-OsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = 'acrolein(T)(28460)',
    structure = SMILES('C=C[CH][O]'),
    E0 = (175.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,338.192,338.217],'cm^-1')),
        HinderedRotor(inertia=(0.430101,'amu*angstrom^2'), symmetry=1, barrier=(34.909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53857,0.0281494,-1.88419e-05,5.72101e-09,-5.0049e-13,21219.9,13.001], Tmin=(100,'K'), Tmax=(1228.78,'K')), NASAPolynomial(coeffs=[8.45062,0.0129027,-5.11122e-06,9.19807e-10,-6.24772e-14,19465.1,-17.9676], Tmin=(1228.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""acrolein(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([O])[CH][O](1973)',
    structure = SMILES('[CH2]C([O])[CH][O]'),
    E0 = (395.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,1187.5,1189.03],'cm^-1')),
        HinderedRotor(inertia=(0.246357,'amu*angstrom^2'), symmetry=1, barrier=(5.66422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24653,'amu*angstrom^2'), symmetry=1, barrier=(5.6682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7725,0.0562038,-9.65776e-05,8.93099e-08,-3.15186e-11,47596.6,21.2675], Tmin=(100,'K'), Tmax=(864.027,'K')), NASAPolynomial(coeffs=[5.58221,0.0227358,-1.09917e-05,2.06761e-09,-1.39916e-13,47529.2,6.86426], Tmin=(864.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][CH][O](28461)',
    structure = SMILES('[CH2][CH][CH][O]'),
    E0 = (538.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1420.43,1420.44],'cm^-1')),
        HinderedRotor(inertia=(0.00413077,'amu*angstrom^2'), symmetry=1, barrier=(5.91441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00413183,'amu*angstrom^2'), symmetry=1, barrier=(5.9158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51471,0.0189427,-7.60192e-06,8.89651e-10,3.10239e-14,64748.2,17.3588], Tmin=(100,'K'), Tmax=(2333.45,'K')), NASAPolynomial(coeffs=[13.7312,0.006622,-3.01959e-06,5.34077e-10,-3.3047e-14,58566.7,-43.6156], Tmin=(2333.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCJCO) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = 'HC2(2881)',
    structure = SMILES('[C]#C'),
    E0 = (556.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0139907,-3.08144e-05,3.10836e-08,-1.10946e-11,66983,5.75946], Tmin=(100,'K'), Tmax=(918.723,'K')), NASAPolynomial(coeffs=[3.14385,0.00498487,-2.32624e-06,4.08757e-10,-2.55787e-14,67315.6,7.08554], Tmin=(918.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.809,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""HC2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([CH][O])O[O](28462)',
    structure = SMILES('[CH2]C([CH][O])O[O]'),
    E0 = (390.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,1575.22],'cm^-1')),
        HinderedRotor(inertia=(0.242886,'amu*angstrom^2'), symmetry=1, barrier=(5.58442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243149,'amu*angstrom^2'), symmetry=1, barrier=(5.59048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243452,'amu*angstrom^2'), symmetry=1, barrier=(5.59743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25979,0.0713595,-0.00013278,1.26771e-07,-4.51592e-11,47072.9,26.1818], Tmin=(100,'K'), Tmax=(877.377,'K')), NASAPolynomial(coeffs=[5.18013,0.02757,-1.36074e-05,2.5527e-09,-1.71262e-13,47382.5,13.4659], Tmin=(877.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(ROOJ) + radical(CJCOOH) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=[C]OOC([CH2])=C[O](25805)',
    structure = SMILES('[CH]=[C]OOC([CH2])=C[O]'),
    E0 = (599.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,890.454],'cm^-1')),
        HinderedRotor(inertia=(0.617143,'amu*angstrom^2'), symmetry=1, barrier=(14.1893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108586,'amu*angstrom^2'), symmetry=1, barrier=(14.2584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616483,'amu*angstrom^2'), symmetry=1, barrier=(14.1742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23766,'amu*angstrom^2'), symmetry=1, barrier=(28.4562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.423368,0.0769347,-9.28101e-05,5.54561e-08,-1.28901e-11,72267.3,34.3643], Tmin=(100,'K'), Tmax=(1058.48,'K')), NASAPolynomial(coeffs=[16.4956,0.0161982,-6.73928e-06,1.24613e-09,-8.64529e-14,68864.9,-44.0884], Tmin=(1058.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(O)CJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]#COOC([CH2])[CH][O](28463)',
    structure = SMILES('[C]#COOC([CH2])[CH][O]'),
    E0 = (912.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2175,525,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.151387,0.105949,-0.000197653,1.84415e-07,-6.4607e-11,109896,33.6056], Tmin=(100,'K'), Tmax=(871.008,'K')), NASAPolynomial(coeffs=[8.3407,0.0341568,-1.7541e-05,3.33253e-09,-2.25001e-13,109661,0.949882], Tmin=(871.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(912.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CJCOOH) + radical(CCsJOH) + radical(Acetyl)"""),
)

species(
    label = 'C#COO[C](C)[CH][O](28464)',
    structure = SMILES('C#COO[C](C)[CH][O]'),
    E0 = (547.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0167035,0.103617,-0.000192142,1.86656e-07,-6.84364e-11,65979.1,30.1909], Tmin=(100,'K'), Tmax=(850.114,'K')), NASAPolynomial(coeffs=[5.21564,0.0434128,-2.28471e-05,4.43698e-09,-3.05467e-13,66386.7,13.5495], Tmin=(850.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(C2CsJOOC) + radical(CCsJOH)"""),
)

species(
    label = 'C#COO[C]([CH2])C[O](28465)',
    structure = SMILES('C#COO[C]([CH2])C[O]'),
    E0 = (581.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180,1317.6],'cm^-1')),
        HinderedRotor(inertia=(0.16828,'amu*angstrom^2'), symmetry=1, barrier=(3.86909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18212,'amu*angstrom^2'), symmetry=1, barrier=(50.1713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166959,'amu*angstrom^2'), symmetry=1, barrier=(3.83871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18282,'amu*angstrom^2'), symmetry=1, barrier=(50.1873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18321,'amu*angstrom^2'), symmetry=1, barrier=(50.1964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184203,0.0980349,-0.000175677,1.69102e-07,-6.2116e-11,70023.8,32.0354], Tmin=(100,'K'), Tmax=(840.704,'K')), NASAPolynomial(coeffs=[5.44766,0.0425217,-2.2264e-05,4.33673e-09,-2.99963e-13,70215.6,13.9596], Tmin=(840.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(C2CsJOOC) + radical(CJCOOH)"""),
)

species(
    label = 'C#COO[C]([CH2])[CH]O(28466)',
    structure = SMILES('C#COO[C]([CH2])[CH]O'),
    E0 = (535.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.406926,0.109004,-0.000192539,1.73169e-07,-5.95769e-11,64585.9,32.8479], Tmin=(100,'K'), Tmax=(854.141,'K')), NASAPolynomial(coeffs=[10.6996,0.0330083,-1.69596e-05,3.25314e-09,-2.22125e-13,63563.4,-13.8624], Tmin=(854.141,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C2CsJOOC) + radical(CCsJOH) + radical(CJCOOH)"""),
)

species(
    label = '[C]#COOC(C)[CH][O](28467)',
    structure = SMILES('[C]#COOC(C)[CH][O]'),
    E0 = (698.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2175,525,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0659325,0.103001,-0.000183753,1.71019e-07,-6.0469e-11,84160.4,31.3034], Tmin=(100,'K'), Tmax=(860.378,'K')), NASAPolynomial(coeffs=[7.52137,0.0381916,-1.92717e-05,3.67174e-09,-2.49775e-13,83948,2.19285], Tmin=(860.378,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(698.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CCsJOH) + radical(Acetyl)"""),
)

species(
    label = '[C]#COOC([CH2])C[O](28468)',
    structure = SMILES('[C]#COOC([CH2])C[O]'),
    E0 = (732.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102168,0.0974078,-0.00016723,1.53348e-07,-5.40733e-11,88205.1,33.1459], Tmin=(100,'K'), Tmax=(850.355,'K')), NASAPolynomial(coeffs=[7.76173,0.0372861,-1.86802e-05,3.56948e-09,-2.44104e-13,87773.5,2.55625], Tmin=(850.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCOJ) + radical(CJCOOH) + radical(Acetyl)"""),
)

species(
    label = '[C]#COOC([CH2])[CH]O(28469)',
    structure = SMILES('[C]#COOC([CH2])[CH]O'),
    E0 = (686.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.495975,0.108473,-0.000184494,1.58063e-07,-5.18829e-11,82767.5,33.9829], Tmin=(100,'K'), Tmax=(871.003,'K')), NASAPolynomial(coeffs=[13.0112,0.0277761,-1.33774e-05,2.48622e-09,-1.66289e-13,81122.6,-25.2516], Tmin=(871.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(686.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsCt) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOH) + radical(CJCOOH) + radical(Acetyl)"""),
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
    E0 = (586.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (575.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (916.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1036.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1067.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1021.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (580.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (639.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (663.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (598.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (598.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (602.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (599.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (706.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (604.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (600.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (752.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (575.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (725.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (644.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (575.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (578.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (879.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (947.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (811.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1124.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (716.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (739.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (685.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (770.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (794.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (766.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['ketene(1375)', 'C=C([O])[CH][O](2850)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.29914e+17,'s^-1'), n=-1.73308, Ea=(10.9298,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.16692063000437535, var=8.814681701569807, Tref=1000.0, N=65, correlation='Root',), comment="""BM rule fitted to 2 training reactions at node Root
    Total Standard Deviation in ln(k): 6.371362717619001
Exact match found for rate rule [Root]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['C#COOC([CH2])C=O(28447)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(20)', '[CH]=[C]OOC=C[O](11850)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][O](1548)', '[CH]=[C]OOC=C(2667)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#COOC([CH2])[C][O](28448)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C([CH][O])OOC#C(28449)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['C#COOC1CC1[O](28450)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.01367e+13,'s^-1'), n=0.0238333, Ea=(5.50893,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_2H;Ypri_rad_out] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['HCCO(2227)', '[O][CH]C1CO1(28451)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.98e+12,'s^-1','*|/',1.2), n=0, Ea=(63.953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for R2OO_S;C_pri_rad_intra;OOR
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['HCCO(2227)', '[CH2]C1OC1[O](23880)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.57916e+14,'s^-1'), n=-0.33125, Ea=(87.6317,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['C#COOC(C)=C[O](28452)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['C#COOC(=C)C[O](28453)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['[O][CH]C1CC=[C]OO1(28454)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.44422e+06,'s^-1'), n=1.00668, Ea=(26.8662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SSS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['[CH2]C1OO[C]=CC1[O](28455)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(315087,'s^-1'), n=1.18716, Ea=(24.4501,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra_cs] for rate rule [R6_SSS_T;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['[CH2]C1[CH]OC=[C]OO1(28456)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['[CH]=C1CC([CH][O])OO1(28457)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.11623e+06,'s^-1'), n=1.18545, Ea=(29.2713,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SSS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['[CH]=C1OOC([CH2])C1[O](28458)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(110269,'s^-1'), n=1.50649, Ea=(24.9081,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra_cs] for rate rule [R6_SSS_T;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['[CH]=C1O[CH]C([CH2])OO1(28459)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.05e+09,'s^-1'), n=0.155, Ea=(177.192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;triplebond_intra_H;radadd_intra] for rate rule [R7;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', 'C#COOC(=C)[CH][O](24131)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(19.2102,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 16.0 to 19.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH][O](1548)', 'C#COOC=C(2350)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(53.1363,'m^3/(mol*s)'), n=1.6135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;Y_1centerbirad] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -4.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(T)(20)', 'C#COOC=C[O](11169)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#CO[O](2880)', 'acrolein(T)(28460)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.91057e-05,'m^3/(mol*s)'), n=3.20111, Ea=(58.1196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;OJ-O2s]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 55.5 to 58.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCCO(2227)', '[CH2]C([O])[CH][O](1973)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(18.3521,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CO[O](2880)', '[CH2][CH][CH][O](28461)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HC2(2881)', '[CH2]C([CH][O])O[O](28462)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]=[C]OOC([CH2])=C[O](25805)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[C]#COOC([CH2])[CH][O](28463)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['C#COO[C](C)[CH][O](28464)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#COO[C]([CH2])C[O](28465)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#COOC([CH2])[CH][O](25545)'],
    products = ['C#COO[C]([CH2])[CH]O(28466)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.23689e+09,'s^-1'), n=1.09705, Ea=(110.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;XH_out] for rate rule [R3HJ;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[C]#COOC(C)[CH][O](28467)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.17493e+06,'s^-1'), n=1.01509, Ea=(71.5445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_TSSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[C]#COOC([CH2])C[O](28468)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(26171.5,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_TSSSS;Ct_rad_out;XH_out]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]#COOC([CH2])[CH]O(28469)'],
    products = ['C#COOC([CH2])[CH][O](25545)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(364667,'s^-1'), n=1.22214, Ea=(79.2357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;Y_rad_out;XH_out] for rate rule [R7HJ_5;Ct_rad_out;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #6189',
    isomers = [
        'C#COOC([CH2])[CH][O](25545)',
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
    label = 'PDepNetwork #6189',
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


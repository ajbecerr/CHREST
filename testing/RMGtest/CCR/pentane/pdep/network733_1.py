species(
    label = 'C#COOC[CH2](2331)',
    structure = SMILES('C#COOC[CH2]'),
    E0 = (365.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.69577,'amu*angstrom^2'), symmetry=1, barrier=(38.989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0558,'amu*angstrom^2'), symmetry=1, barrier=(24.2749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69092,'amu*angstrom^2'), symmetry=1, barrier=(38.8777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05729,'amu*angstrom^2'), symmetry=1, barrier=(24.3093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43382,0.0599935,-6.79368e-05,4.16907e-08,-1.05232e-11,44073.8,22.1748], Tmin=(100,'K'), Tmax=(950.206,'K')), NASAPolynomial(coeffs=[9.87748,0.0244488,-1.18257e-05,2.32296e-09,-1.65503e-13,42469.2,-18.1297], Tmin=(950.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCOOH)"""),
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
    label = 'C#COO[CH2](2887)',
    structure = SMILES('C#COO[CH2]'),
    E0 = (382.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.81868,'amu*angstrom^2'), symmetry=1, barrier=(41.8149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81109,'amu*angstrom^2'), symmetry=1, barrier=(41.6406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22571,0.041128,-4.32669e-05,2.39862e-08,-5.46916e-12,46113.1,15.1973], Tmin=(100,'K'), Tmax=(1043.81,'K')), NASAPolynomial(coeffs=[8.57818,0.0167845,-8.28394e-06,1.64292e-09,-1.17754e-13,44786.9,-15.722], Tmin=(1043.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-OsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CsJOOC)"""),
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
    label = '[CH]COOC#C(2888)',
    structure = SMILES('[CH]COOC#C'),
    E0 = (599.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100,180,451,1356.83],'cm^-1')),
        HinderedRotor(inertia=(0.0684747,'amu*angstrom^2'), symmetry=1, barrier=(17.4219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79985,'amu*angstrom^2'), symmetry=1, barrier=(41.3821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79962,'amu*angstrom^2'), symmetry=1, barrier=(41.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4279,0.0605373,-7.49608e-05,4.96004e-08,-1.33548e-11,72247.7,21.2059], Tmin=(100,'K'), Tmax=(897.883,'K')), NASAPolynomial(coeffs=[10.0044,0.0223282,-1.11266e-05,2.20268e-09,-1.57282e-13,70707.6,-19.2468], Tmin=(897.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCJ2_triplet)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'oxirane(2889)',
    structure = SMILES('C1CO1'),
    E0 = (-62.2743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85636,-0.00852312,7.25147e-05,-8.57001e-08,3.15182e-11,-7474.05,7.20137], Tmin=(100,'K'), Tmax=(947.618,'K')), NASAPolynomial(coeffs=[6.87395,0.00877628,-2.41484e-06,4.63775e-10,-3.81439e-14,-9394.58,-14.3103], Tmin=(947.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.2743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""oxirane""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[C]1=CCCOO1(2890)',
    structure = SMILES('[C]1=CCCOO1'),
    E0 = (201.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54467,0.0223409,2.86489e-05,-5.1602e-08,2.1424e-11,24351.8,17.3932], Tmin=(100,'K'), Tmax=(936.269,'K')), NASAPolynomial(coeffs=[8.89938,0.0196209,-6.13142e-06,1.03118e-09,-7.12144e-14,22091.1,-18.5642], Tmin=(936.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1CCOO1(2891)',
    structure = SMILES('[CH]=C1CCOO1'),
    E0 = (226.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55865,0.0160364,5.60877e-05,-8.32054e-08,3.23746e-11,27350.7,19.8054], Tmin=(100,'K'), Tmax=(971.023,'K')), NASAPolynomial(coeffs=[11.8733,0.0169837,-6.11234e-06,1.1982e-09,-9.19218e-14,23688.1,-34.4032], Tmin=(971.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Cds_P)"""),
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
    label = 'C2H4(22)',
    structure = SMILES('C=C'),
    E0 = (41.9072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98876,-0.00674712,5.04398e-05,-5.70744e-08,2.04945e-11,5047.05,3.80494], Tmin=(100,'K'), Tmax=(946.013,'K')), NASAPolynomial(coeffs=[4.59022,0.00872722,-2.66493e-06,4.8171e-10,-3.6069e-14,4127.02,-3.3246], Tmin=(946.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.9072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[O](1358)',
    structure = SMILES('[CH2]C[O]'),
    E0 = (187.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1795.79],'cm^-1')),
        HinderedRotor(inertia=(0.00405474,'amu*angstrom^2'), symmetry=1, barrier=(9.27899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63182,0.0106853,2.91824e-06,-4.98315e-09,1.21781e-12,22581.3,11.2346], Tmin=(100,'K'), Tmax=(1659.64,'K')), NASAPolynomial(coeffs=[5.14351,0.0138609,-6.11494e-06,1.12108e-09,-7.46163e-14,21140.5,0.346415], Tmin=(1659.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = 'C2H4(T)(23)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,995.342,995.492,2246.13,2948.67],'cm^-1')),
        HinderedRotor(inertia=(0.00723049,'amu*angstrom^2'), symmetry=1, barrier=(5.08405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40735,0.0100313,6.40897e-06,-1.41287e-08,5.92653e-12,38288.2,6.11706], Tmin=(100,'K'), Tmax=(954.266,'K')), NASAPolynomial(coeffs=[5.52251,0.00856169,-2.90741e-06,5.02347e-10,-3.44567e-14,37547.7,-5.75289], Tmin=(954.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]CO[O](49)',
    structure = SMILES('[CH2]CO[O]'),
    E0 = (188.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.345977,'amu*angstrom^2'), symmetry=1, barrier=(7.9547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000700344,'amu*angstrom^2'), symmetry=1, barrier=(7.95174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81264,0.0259367,-2.20523e-05,1.10177e-08,-2.31287e-12,22700.9,17.2501], Tmin=(100,'K'), Tmax=(1124.17,'K')), NASAPolynomial(coeffs=[6.40981,0.0131371,-4.97349e-06,8.8932e-10,-6.04387e-14,21892.2,-0.525236], Tmin=(1124.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CJCOOH)"""),
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
    label = '[C]#COOC[CH2](2878)',
    structure = SMILES('[C]#COOC[CH2]'),
    E0 = (702.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2175,525,223.324],'cm^-1')),
        HinderedRotor(inertia=(1.39019,'amu*angstrom^2'), symmetry=1, barrier=(47.6808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39897,'amu*angstrom^2'), symmetry=1, barrier=(47.7038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263076,'amu*angstrom^2'), symmetry=1, barrier=(8.69415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.58146,'amu*angstrom^2'), symmetry=1, barrier=(89.8734,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19365,0.0668206,-0.000102903,8.60762e-08,-2.84874e-11,84629.6,23.7669], Tmin=(100,'K'), Tmax=(830.922,'K')), NASAPolynomial(coeffs=[8.85828,0.0228854,-1.08843e-05,2.0536e-09,-1.40346e-13,83598.8,-10.3289], Tmin=(830.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(702.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CJCOOH) + radical(Acetyl)"""),
)

species(
    label = 'C#COO[CH]C(2892)',
    structure = SMILES('C#COO[CH]C'),
    E0 = (338.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.54039,'amu*angstrom^2'), symmetry=1, barrier=(35.4165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54024,'amu*angstrom^2'), symmetry=1, barrier=(35.4131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54014,'amu*angstrom^2'), symmetry=1, barrier=(35.4108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5403,'amu*angstrom^2'), symmetry=1, barrier=(35.4145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53523,0.0587323,-6.73272e-05,4.39914e-08,-1.2096e-11,40759.1,20.1142], Tmin=(100,'K'), Tmax=(867.742,'K')), NASAPolynomial(coeffs=[8.28483,0.027618,-1.35407e-05,2.66727e-09,-1.90033e-13,39587.8,-11.491], Tmin=(867.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOOC)"""),
)

species(
    label = '[C]#COOCC(2893)',
    structure = SMILES('[C]#COOCC'),
    E0 = (488.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,180],'cm^-1')),
        HinderedRotor(inertia=(2.72936,'amu*angstrom^2'), symmetry=1, barrier=(62.7534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.556092,'amu*angstrom^2'), symmetry=1, barrier=(12.7856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28734,'amu*angstrom^2'), symmetry=1, barrier=(29.5986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72676,'amu*angstrom^2'), symmetry=1, barrier=(62.6937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31707,0.0633892,-8.71527e-05,7.00661e-08,-2.31518e-11,58892.1,21.3308], Tmin=(100,'K'), Tmax=(778.278,'K')), NASAPolynomial(coeffs=[7.90524,0.0271591,-1.27578e-05,2.42742e-09,-1.68047e-13,57938.4,-8.34075], Tmin=(778.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Acetyl)"""),
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
    E0 = (404.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (820.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (811.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (429.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (392.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (394.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (513.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (429.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (365.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (659.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (745.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (761.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (914.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (526.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (560.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#COOC[CH2](2331)'],
    products = ['ketene(1375)', 'vinoxy(1351)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.59828e+17,'s^-1'), n=-1.73308, Ea=(39.1142,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.16692063000437535, var=8.814681701569807, Tref=1000.0, N=65, correlation='Root',), comment="""BM rule fitted to 2 training reactions at node Root
    Total Standard Deviation in ln(k): 6.371362717619001
Exact match found for rate rule [Root]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Retroene"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(20)', 'C#COO[CH2](2887)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH]COOC#C(2888)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#COOC[CH2](2331)'],
    products = ['HCCO(2227)', 'oxirane(2889)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.98e+12,'s^-1','*|/',1.2), n=0, Ea=(63.953,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for R2OO_S;C_pri_rad_intra;OOR
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#COOC[CH2](2331)'],
    products = ['[C]1=CCCOO1(2890)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.44422e+06,'s^-1'), n=1.00668, Ea=(26.8662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SSS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C#COOC[CH2](2331)'],
    products = ['[CH]=C1CCOO1(2891)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.11623e+06,'s^-1'), n=1.18545, Ea=(29.2713,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SSS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C#COOC=C(2350)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#CO[O](2880)', 'C2H4(22)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(71.2,'cm^3/(mol*s)'), n=3.22, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), Tmin=(400,'K'), Tmax=(1100,'K'), comment="""From training reaction 2769 used for Cds-HH_Cds-HH;OJ-O2s
Exact match found for rate rule [Cds-HH_Cds-HH;OJ-O2s]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['HCCO(2227)', '[CH2]C[O](1358)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(12.8556,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
family: R_Recombination
Ea raised from 12.7 to 12.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CO[O](2880)', 'C2H4(T)(23)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(212954,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['HC2(2881)', '[CH2]CO[O](49)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.36745e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]=[C]OOC=C(2667)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.82867e+07,'m^3/(mol*s)'), n=0.0631113, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0175378549852, var=0.221368827459, Tref=1000.0, N=8, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN
    Total Standard Deviation in ln(k): 0.987289785558
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[C]#COOC[CH2](2878)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.82867e+07,'m^3/(mol*s)'), n=0.0631113, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0175378549852, var=0.221368827459, Tref=1000.0, N=8, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN
    Total Standard Deviation in ln(k): 0.987289785558
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#COO[CH]C(2892)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.41e+13,'s^-1','+|-',2), n=0, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1500,'K'), comment="""From training reaction 346 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#COOCC(2893)'],
    products = ['C#COOC[CH2](2331)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.17493e+06,'s^-1'), n=1.01509, Ea=(71.5445,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_TSSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #733',
    isomers = [
        'C#COOC[CH2](2331)',
    ],
    reactants = [
        ('ketene(1375)', 'vinoxy(1351)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #733',
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


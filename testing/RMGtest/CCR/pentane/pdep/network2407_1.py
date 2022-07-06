species(
    label = 'C=C([O])[CH][C]=O(2676)',
    structure = SMILES('[CH2]C([O])=C[C]=O'),
    E0 = (125.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.67906,'amu*angstrom^2'), symmetry=1, barrier=(15.6129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675636,'amu*angstrom^2'), symmetry=1, barrier=(15.5342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36779,0.0612201,-9.15025e-05,6.93676e-08,-2.06256e-11,15176.3,21.4686], Tmin=(100,'K'), Tmax=(828.259,'K')), NASAPolynomial(coeffs=[11.0334,0.0145441,-6.97656e-06,1.33711e-09,-9.28003e-14,13575,-23.3417], Tmin=(828.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CCJ=O)"""),
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
    label = '[CH2][C]=C[C]=O(9897)',
    structure = SMILES('[CH2][C]=C[C]=O'),
    E0 = (432.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.00399,'amu*angstrom^2'), symmetry=1, barrier=(23.0837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00442,'amu*angstrom^2'), symmetry=1, barrier=(23.0936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33986,0.0429161,-4.87447e-05,1.23757e-08,1.41321e-11,52082.9,16.5479], Tmin=(100,'K'), Tmax=(512.282,'K')), NASAPolynomial(coeffs=[6.6704,0.019476,-1.04855e-05,2.11536e-09,-1.51263e-13,51503.1,-2.77645], Tmin=(512.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Cds_S) + radical(C=CCJ=O)"""),
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
    label = '[O][C]=C[C]=O(9898)',
    structure = SMILES('O=[C][CH][C]=O'),
    E0 = (211.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1850,1860,440,470,900,1000],'cm^-1')),
        HinderedRotor(inertia=(1.2407,'amu*angstrom^2'), symmetry=1, barrier=(28.5261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24145,'amu*angstrom^2'), symmetry=1, barrier=(28.5433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0388,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81165,0.0557507,-0.00011026,1.00914e-07,-3.40854e-11,25521.6,14.8862], Tmin=(100,'K'), Tmax=(897.886,'K')), NASAPolynomial(coeffs=[7.74751,0.0112957,-5.90386e-06,1.09028e-09,-7.08098e-14,25181.6,-9.06853], Tmin=(897.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCJC=O) + radical(C=OCCJ=O) + radical(C=OCCJ=O)"""),
)

species(
    label = '[C]=O(2355)',
    structure = SMILES('[C]=O'),
    E0 = (438.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3323.79],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09152,0.00193295,-1.59153e-05,2.47563e-08,-1.11287e-11,52770.5,4.46624], Tmin=(100,'K'), Tmax=(866.029,'K')), NASAPolynomial(coeffs=[1.05092,0.00526657,-3.13864e-06,6.40624e-10,-4.47939e-14,53698.8,21.0169], Tmin=(866.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
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
    label = '[CH]=C([O])C=[C][O](9624)',
    structure = SMILES('[CH]C([O])=C[C]=O'),
    E0 = (337.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09785,'amu*angstrom^2'), symmetry=1, barrier=(48.2337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10192,'amu*angstrom^2'), symmetry=1, barrier=(48.3274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46965,0.0598019,-8.50278e-05,6.58157e-08,-2.06105e-11,40641.6,21.5237], Tmin=(100,'K'), Tmax=(778.943,'K')), NASAPolynomial(coeffs=[9.00533,0.0211058,-1.05129e-05,2.04269e-09,-1.43112e-13,39467.6,-12.9492], Tmin=(778.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=[C]C=C1CO1(9899)',
    structure = SMILES('O=[C]C=C1CO1'),
    E0 = (61.4183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56162,0.0455057,-3.58585e-05,6.31676e-09,2.52131e-12,7481.88,16.8268], Tmin=(100,'K'), Tmax=(1010.34,'K')), NASAPolynomial(coeffs=[15.5189,0.00670776,-2.69421e-06,5.58114e-10,-4.36296e-14,3821.44,-54.8103], Tmin=(1010.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.4183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(methyleneoxirane) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C1=CC(=O)O1(9900)',
    structure = SMILES('C=C1[CH]C(=O)O1'),
    E0 = (-127.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11752,0.0230823,4.49887e-05,-8.49632e-08,3.68134e-11,-15295.1,15.9222], Tmin=(100,'K'), Tmax=(937.084,'K')), NASAPolynomial(coeffs=[17.4164,0.00440112,2.61788e-07,-4.93883e-11,-4.95417e-15,-20209.4,-67.8146], Tmin=(937.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(4-Methylene-2-oxetanone) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C1=CC(=O)C1(9901)',
    structure = SMILES('[O]C1=CC(=O)C1'),
    E0 = (-37.7124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63788,0.0179248,3.07636e-05,-5.00229e-08,1.86844e-11,-4476.34,15.0334], Tmin=(100,'K'), Tmax=(1049.41,'K')), NASAPolynomial(coeffs=[11.3583,0.0147358,-7.6314e-06,1.65605e-09,-1.27771e-13,-7961.24,-35.3418], Tmin=(1049.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-37.7124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=C=C=O(9902)',
    structure = SMILES('C=C(O)[C]=C=O'),
    E0 = (44.5108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33834,0.0564656,-7.333e-05,4.65222e-08,-1.13064e-11,5451.07,19.1039], Tmin=(100,'K'), Tmax=(1021.81,'K')), NASAPolynomial(coeffs=[13.6363,0.00832452,-2.66054e-06,4.15331e-10,-2.58614e-14,2937.82,-40.492], Tmin=(1021.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.5108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(Cds_S)"""),
)

species(
    label = 'CC([O])=C=C=O(9903)',
    structure = SMILES('CC(=O)[C]=C=O'),
    E0 = (36.2561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34652,0.0405903,-4.86098e-05,3.73927e-08,-1.28457e-11,4416.25,19.4257], Tmin=(100,'K'), Tmax=(687.37,'K')), NASAPolynomial(coeffs=[5.17234,0.0241462,-1.27253e-05,2.58943e-09,-1.87721e-13,4027.77,6.85201], Tmin=(687.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.2561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=CJC=O)"""),
)

species(
    label = '[CH2][C]1OC1[C]=O(9904)',
    structure = SMILES('[CH2][C]1OC1[C]=O'),
    E0 = (331.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27179,0.04711,-5.10079e-05,2.7904e-08,-5.5465e-12,39990.4,24.8002], Tmin=(100,'K'), Tmax=(1522.36,'K')), NASAPolynomial(coeffs=[11.4275,0.00737047,1.01156e-06,-5.09306e-10,4.45565e-14,38411.2,-23.4946], Tmin=(1522.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1CC1[C]=O(9905)',
    structure = SMILES('[O][C]1CC1[C]=O'),
    E0 = (322.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03172,0.0376823,-2.66889e-05,5.35198e-09,1.21113e-12,38916.4,20.864], Tmin=(100,'K'), Tmax=(1031.22,'K')), NASAPolynomial(coeffs=[11.2701,0.012201,-4.68435e-06,8.62591e-10,-6.08625e-14,36460.5,-26.659], Tmin=(1031.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1([O])[CH]C1=O(9906)',
    structure = SMILES('[CH2]C1([O])C=C1[O]'),
    E0 = (409.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47184,0.0489822,-4.49124e-05,1.52603e-08,-3.02672e-13,49375.7,19.5962], Tmin=(100,'K'), Tmax=(988.775,'K')), NASAPolynomial(coeffs=[14.9837,0.00772894,-2.66975e-06,4.92562e-10,-3.61587e-14,46048.3,-48.7522], Tmin=(988.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)2OJ) + radical(C=C(C)OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'CO(2039)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(815.652,'J/mol'), sigma=(3.65,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.95,'angstroms^3'), rotrelaxcollnum=1.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C([O])=C=C=O(9907)',
    structure = SMILES('C=C([O])[C]=C=O'),
    E0 = (182.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.58125,'amu*angstrom^2'), symmetry=1, barrier=(13.3641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74143,0.0520469,-7.92674e-05,6.15818e-08,-1.8499e-11,22006.6,19.1592], Tmin=(100,'K'), Tmax=(896.992,'K')), NASAPolynomial(coeffs=[9.75664,0.011928,-4.86006e-06,8.41339e-10,-5.41586e-14,20744.7,-17.6569], Tmin=(896.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]C=C=O(9594)',
    structure = SMILES('O=[C]C=C=O'),
    E0 = (-7.87925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1855,455,950,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.841576,'amu*angstrom^2'), symmetry=1, barrier=(19.3495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0388,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98472,0.0259268,-2.97918e-05,2.05261e-08,-6.41245e-12,-914.289,13.446], Tmin=(100,'K'), Tmax=(741.746,'K')), NASAPolynomial(coeffs=[4.88246,0.0156926,-9.09529e-06,1.92406e-09,-1.42625e-13,-1195.81,4.85747], Tmin=(741.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.87925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + missing(Cdd-CdO2d) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C=C[C]=O(9908)',
    structure = SMILES('C=C=C[C]=O'),
    E0 = (213.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.13978,'amu*angstrom^2'), symmetry=1, barrier=(26.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00684,0.0284132,-1.9291e-05,6.01948e-09,-7.90367e-13,25757,13.5223], Tmin=(100,'K'), Tmax=(1551.73,'K')), NASAPolynomial(coeffs=[6.76337,0.0187298,-9.93048e-06,1.99793e-09,-1.42457e-13,24591.2,-6.25135], Tmin=(1551.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]C([O])=[C][C]=O(9909)',
    structure = SMILES('[CH2]C([O])=[C][C]=O'),
    E0 = (369.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.680006,'amu*angstrom^2'), symmetry=1, barrier=(15.6347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680123,'amu*angstrom^2'), symmetry=1, barrier=(15.6374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15508,0.0670571,-0.000114989,9.48772e-08,-3.01365e-11,44524.5,22.9321], Tmin=(100,'K'), Tmax=(829.692,'K')), NASAPolynomial(coeffs=[11.9656,0.0108003,-5.80085e-06,1.13136e-09,-7.79261e-14,42873,-26.3458], Tmin=(829.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C([O])=[C]C=O(9910)',
    structure = SMILES('C=C([O])[C]=C[O]'),
    E0 = (149.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55298,'amu*angstrom^2'), symmetry=1, barrier=(35.706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562656,0.0598018,-6.95646e-05,3.74263e-08,-7.30439e-12,18117.9,21.9387], Tmin=(100,'K'), Tmax=(1502.94,'K')), NASAPolynomial(coeffs=[17.396,0.00121988,2.6563e-06,-7.09798e-10,5.39892e-14,14614.4,-60.9531], Tmin=(1502.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(O)=[C][C]=O(9911)',
    structure = SMILES('[CH2]C(O)=[C][C]=O'),
    E0 = (231.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.836731,'amu*angstrom^2'), symmetry=1, barrier=(19.2381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837732,'amu*angstrom^2'), symmetry=1, barrier=(19.2611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83643,'amu*angstrom^2'), symmetry=1, barrier=(19.2312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897348,0.0696606,-0.000102198,7.01207e-08,-1.84032e-11,27962.8,22.3619], Tmin=(100,'K'), Tmax=(944.772,'K')), NASAPolynomial(coeffs=[15.515,0.0077724,-3.93902e-06,7.86111e-10,-5.63883e-14,25200.7,-47.3293], Tmin=(944.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(O)CJ) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC([O])=[C][C]=O(9912)',
    structure = SMILES('CC([O])=[C][C]=O'),
    E0 = (210.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.650692,'amu*angstrom^2'), symmetry=1, barrier=(14.9607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652899,'amu*angstrom^2'), symmetry=1, barrier=(15.0114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24499,0.0658984,-0.000109799,9.35559e-08,-3.13656e-11,25407.2,21.9193], Tmin=(100,'K'), Tmax=(782.833,'K')), NASAPolynomial(coeffs=[10.0482,0.0171777,-9.27892e-06,1.8503e-09,-1.30433e-13,24143.5,-17.6638], Tmin=(782.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
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
    label = 'O=[C]C1CC1=O(9877)',
    structure = SMILES('O=[C]C1CC1=O'),
    E0 = (48.4459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5683,0.0198976,2.80401e-05,-5.33332e-08,2.2539e-11,5888.79,18.5875], Tmin=(100,'K'), Tmax=(951.646,'K')), NASAPolynomial(coeffs=[11.6685,0.0106739,-3.17356e-06,5.84395e-10,-4.52484e-14,2842.4,-31.7703], Tmin=(951.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.4459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C][CH][C]1CO1(9913)',
    structure = SMILES('O=[C][CH][C]1CO1'),
    E0 = (329.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92601,0.0466042,-6.33752e-05,4.68236e-08,-1.28603e-11,39758.6,22.9002], Tmin=(100,'K'), Tmax=(1102.04,'K')), NASAPolynomial(coeffs=[7.39791,0.0144972,-3.00572e-06,2.20432e-10,-8.42282e-16,39296.1,-0.656326], Tmin=(1102.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]1[CH]C(=O)O1(9914)',
    structure = SMILES('[CH2][C]1[CH]C(=O)O1'),
    E0 = (277.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4601,0.0279092,-2.22959e-07,-2.37035e-08,1.31001e-11,33465.5,22.6806], Tmin=(100,'K'), Tmax=(879.897,'K')), NASAPolynomial(coeffs=[9.67505,0.0122473,-2.73771e-06,3.36631e-10,-1.92789e-14,31532.4,-14.9739], Tmin=(879.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C2CsJOC(O)) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C1OC1[C]=O(9885)',
    structure = SMILES('C=C1OC1[C]=O'),
    E0 = (60.9201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03026,0.0246384,3.72491e-05,-8.18572e-08,3.75728e-11,7415.29,17.3107], Tmin=(100,'K'), Tmax=(920.955,'K')), NASAPolynomial(coeffs=[19.7401,-0.00374361,4.42112e-06,-8.66589e-10,5.26248e-14,2094.93,-77.8459], Tmin=(920.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.9201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1[CH]C(=O)C1(9915)',
    structure = SMILES('[O][C]1C=C([O])C1'),
    E0 = (301.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07154,0.0306602,6.65968e-06,-3.89868e-08,1.94486e-11,36389.9,19.5866], Tmin=(100,'K'), Tmax=(937.324,'K')), NASAPolynomial(coeffs=[14.9481,0.0054728,-6.622e-07,9.69714e-11,-1.10119e-14,32668.6,-48.6761], Tmin=(937.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH]=C(O)[CH][C]=O(9916)',
    structure = SMILES('[CH]C(O)=C[C]=O'),
    E0 = (199.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.96874,'amu*angstrom^2'), symmetry=1, barrier=(45.2652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97311,'amu*angstrom^2'), symmetry=1, barrier=(45.3657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9666,'amu*angstrom^2'), symmetry=1, barrier=(45.216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17961,0.0628766,-7.4291e-05,4.42506e-08,-1.04445e-11,24081.2,21.0633], Tmin=(100,'K'), Tmax=(1031.96,'K')), NASAPolynomial(coeffs=[12.8955,0.0174644,-8.2824e-06,1.6077e-09,-1.1396e-13,21663.1,-35.8276], Tmin=(1031.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C([O])C[C]=O(2677)',
    structure = SMILES('[CH]=C([O])C[C]=O'),
    E0 = (237.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1855,455,950,3120,650,792.5,1650,457.716],'cm^-1')),
        HinderedRotor(inertia=(0.0811444,'amu*angstrom^2'), symmetry=1, barrier=(11.9485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0819876,'amu*angstrom^2'), symmetry=1, barrier=(12.0127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73346,0.0430516,-3.78896e-05,1.47612e-08,-1.97815e-12,28611.6,22.588], Tmin=(100,'K'), Tmax=(1203.61,'K')), NASAPolynomial(coeffs=[13.8763,0.00944907,-4.42733e-06,8.87773e-10,-6.46617e-14,25199.4,-40.2764], Tmin=(1203.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C=C[O](9605)',
    structure = SMILES('[CH]C([O])=CC=O'),
    E0 = (176.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,2782.5,750,1395,475,1775,1000,263.046,263.054,263.089,263.095],'cm^-1')),
        HinderedRotor(inertia=(1.01488,'amu*angstrom^2'), symmetry=1, barrier=(49.8422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01453,'amu*angstrom^2'), symmetry=1, barrier=(49.842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0356,0.0457448,-3.88501e-05,1.75926e-08,-3.39965e-12,21304.6,19.9636], Tmin=(100,'K'), Tmax=(1182.13,'K')), NASAPolynomial(coeffs=[8.35906,0.024348,-1.16999e-05,2.28112e-09,-1.61542e-13,19809.5,-11.6015], Tmin=(1182.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([O])C=C=O(2675)',
    structure = SMILES('C=C([O])C=C=O'),
    E0 = (-55.5262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.677304,'amu*angstrom^2'), symmetry=1, barrier=(15.5726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76649,0.0475408,-5.47877e-05,3.26419e-08,-7.62371e-12,-6596.45,18.5953], Tmin=(100,'K'), Tmax=(1052.25,'K')), NASAPolynomial(coeffs=[11.0503,0.0122492,-4.47858e-06,7.67665e-10,-5.07754e-14,-8550.21,-26.6665], Tmin=(1052.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.5262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([O])[C]=[C][O](9917)',
    structure = SMILES('[CH2]C([O])[C][C]=O'),
    E0 = (650.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1855,455,950,287.75,297.294,301.116],'cm^-1')),
        HinderedRotor(inertia=(0.0682857,'amu*angstrom^2'), symmetry=1, barrier=(4.12594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0682352,'amu*angstrom^2'), symmetry=1, barrier=(4.17606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382579,'amu*angstrom^2'), symmetry=1, barrier=(23.5144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30892,0.0636388,-0.000102728,8.46844e-08,-2.70261e-11,78336,24.8197], Tmin=(100,'K'), Tmax=(864.52,'K')), NASAPolynomial(coeffs=[10.0061,0.016055,-7.42537e-06,1.36764e-09,-9.14561e-14,77106.6,-14.2858], Tmin=(864.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(650.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1C=[C]OC1(9918)',
    structure = SMILES('[O]C1=C[C]OC1'),
    E0 = (304.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23678,0.0262523,1.77654e-05,-4.88185e-08,2.24465e-11,36715.8,16.8095], Tmin=(100,'K'), Tmax=(941.878,'K')), NASAPolynomial(coeffs=[14.3151,0.00655347,-1.18092e-06,2.07157e-10,-1.94976e-14,33039.1,-48.1777], Tmin=(941.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C1([O])C=[C]O1(9919)',
    structure = SMILES('[CH2]C1([O])C=[C]O1'),
    E0 = (432.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19208,0.0453532,-7.01159e-06,-4.85544e-08,2.9656e-11,52175.1,19.8703], Tmin=(100,'K'), Tmax=(889.415,'K')), NASAPolynomial(coeffs=[23.6896,-0.00974869,8.20845e-06,-1.71521e-09,1.17874e-13,46350.7,-96.2763], Tmin=(889.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CC(C)(O)OJ) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = 'C=C([O])[C]=[C]O(9920)',
    structure = SMILES('C=C([O])[C]=[C]O'),
    E0 = (247.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11612,'amu*angstrom^2'), symmetry=1, barrier=(25.6617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11617,'amu*angstrom^2'), symmetry=1, barrier=(25.663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.712271,0.0645199,-8.64191e-05,5.33794e-08,-1.20339e-11,29926.3,23.1117], Tmin=(100,'K'), Tmax=(1275.32,'K')), NASAPolynomial(coeffs=[17.0096,0.00109766,2.65062e-06,-7.4764e-10,5.97827e-14,26770.2,-55.5535], Tmin=(1275.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])C=[C]O(9921)',
    structure = SMILES('[CH]=C([O])C=[C]O'),
    E0 = (295.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.927454,'amu*angstrom^2'), symmetry=1, barrier=(21.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925422,'amu*angstrom^2'), symmetry=1, barrier=(21.2773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544173,0.0634832,-8.06967e-05,4.62669e-08,-9.6109e-12,35721.5,24.5244], Tmin=(100,'K'), Tmax=(1394.16,'K')), NASAPolynomial(coeffs=[18.6433,-0.00168628,3.66644e-06,-8.86636e-10,6.6248e-14,31961.7,-64.1929], Tmin=(1394.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'O=[C][CH]C[C]=O(2650)',
    structure = SMILES('O=[C][CH]C[C]=O'),
    E0 = (209.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1850,1860,440,470,900,1000,293.378],'cm^-1')),
        HinderedRotor(inertia=(0.121705,'amu*angstrom^2'), symmetry=1, barrier=(7.46704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122394,'amu*angstrom^2'), symmetry=1, barrier=(7.46743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000934063,'amu*angstrom^2'), symmetry=1, barrier=(7.46757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3653.08,'J/mol'), sigma=(5.8998,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.60 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95484,0.0510089,-8.46738e-05,7.67342e-08,-2.66767e-11,25315.5,25.1424], Tmin=(100,'K'), Tmax=(869.639,'K')), NASAPolynomial(coeffs=[5.67682,0.0207036,-9.65797e-06,1.79183e-09,-1.20387e-13,25166.7,10.5726], Tmin=(869.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]1C=[C]OO1(9922)',
    structure = SMILES('C=C1[CH][C]OO1'),
    E0 = (557.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21869,0.0222034,4.06311e-05,-7.4612e-08,3.11944e-11,67151.3,16.5538], Tmin=(100,'K'), Tmax=(964.864,'K')), NASAPolynomial(coeffs=[15.9396,0.00760668,-2.41448e-06,5.51376e-10,-4.94922e-14,62535.2,-59.3509], Tmin=(964.864,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJCO) + radical(CH2_triplet)"""),
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
    label = '[CH2][C]([O])[C]=C[O](9923)',
    structure = SMILES('[CH2]C([O])=[C][CH][O]'),
    E0 = (479.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,605.906,606.205,606.746],'cm^-1')),
        HinderedRotor(inertia=(0.00869473,'amu*angstrom^2'), symmetry=1, barrier=(2.26831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00880927,'amu*angstrom^2'), symmetry=1, barrier=(2.2799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94374,0.0449052,-4.91862e-05,2.88638e-08,-6.78629e-12,57801.1,24.289], Tmin=(100,'K'), Tmax=(1034.11,'K')), NASAPolynomial(coeffs=[9.61839,0.0152189,-6.12504e-06,1.10303e-09,-7.49436e-14,56213.8,-12.9941], Tmin=(1034.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=C(O)CJ) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])[C]=C=O(9924)',
    structure = SMILES('[CH2]C([O])[C]=C=O'),
    E0 = (361.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.686114,'amu*angstrom^2'), symmetry=1, barrier=(15.7751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89243,'amu*angstrom^2'), symmetry=1, barrier=(43.5107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22428,0.0656877,-0.000106325,8.75196e-08,-2.77083e-11,43598,21.2505], Tmin=(100,'K'), Tmax=(885.624,'K')), NASAPolynomial(coeffs=[10.1351,0.0163152,-7.24529e-06,1.30024e-09,-8.52416e-14,42377.5,-18.6363], Tmin=(885.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CC(C)OJ) + radical(CJCO) + radical(CCCJ=C=O)"""),
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
    E0 = (125.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (951.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (649.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (816.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (548.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (130.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (133.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (133.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (188.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (188.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (331.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (322.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (409.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (223.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (410.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (412.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (471.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (581.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (312.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (424.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (366.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (641.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (130.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (329.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (367.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (388.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (128.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (301.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (391.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (429.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (385.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (125.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (673.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (304.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (432.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (348.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (440.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (328.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (455.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (557.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (869.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (502.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (539.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['HCCO(2227)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH2][C]=C[C]=O(9897)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(20)', '[O][C]=C[C]=O(9898)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[C]=O(2355)', '[CH]=C([CH2])[O](9653)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]=C([O])C=[C][O](9624)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['O=[C]C=C1CO1(9899)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2]C1=CC(=O)O1(9900)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[O]C1=CC(=O)C1(9901)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2]C(O)=C=C=O(9902)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['CC([O])=C=C=O(9903)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2][C]1OC1[C]=O(9904)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(206.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 203.2 to 206.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[O][C]1CC1[C]=O(9905)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(197.519,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 194.1 to 197.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2]C1([O])[CH]C1=O(9906)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(284.307,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 282.7 to 284.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['CO(2039)', '[CH]=C([CH2])[O](9653)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2]C([O])=C=C=O(9907)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cdd;HJ] for rate rule [Ca_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2(T)(20)', 'O=[C]C=C=O(9594)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.96577e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;CH2_triplet]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O(4)', 'C=C=C[C]=O(9908)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C([O])=[C][C]=O(9909)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])=[C]C=O(9910)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)=[C][C]=O(9911)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC([O])=[C][C]=O(9912)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C][O](6861)', '[CH2][C]=O(1376)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['O=[C]C1CC1=O(9877)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['O=[C][CH][C]1CO1(9913)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.95361e+10,'s^-1'), n=0.549916, Ea=(204.537,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R3_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 202.1 to 204.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2][C]1[CH]C(=O)O1(9914)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra] for rate rule [R4_linear;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C][O](6861)', 'ketene(1375)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.33566,'m^3/(mol*s)'), n=1.91381, Ea=(23.9557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ck_Cds-HH;YJ] for rate rule [Ck_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['C=C1OC1[C]=O(9885)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[O][C]1[CH]C(=O)C1(9915)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.77275e+09,'s^-1'), n=0.708429, Ea=(176.476,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 172.4 to 176.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(O)[CH][C]=O(9916)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C([O])C[C]=O(2677)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([O])C=C[O](9605)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_2;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['C=C([O])C=C=O(2675)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([O])[C]=[C][O](9917)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[O][C]1C=[C]OC1(9918)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(179.229,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 174.6 to 179.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2]C1([O])C=[C]O1(9919)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(307.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 305.8 to 307.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['HCCO(2227)', '[CH2][C]=O(1376)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C([O])[C]=[C]O(9920)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C([O])C=[C]O(9921)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C][CH]C[C]=O(2650)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([O])[CH][C]=O(2676)'],
    products = ['[CH2][C]1C=[C]OO1(9922)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.44463e+09,'s^-1'), n=0.470283, Ea=(432.247,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 427.4 to 432.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['HCCO(2227)', '[CH2][C][O](2821)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][C]([O])[C]=C[O](9923)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])[C]=C=O(9924)'],
    products = ['C=C([O])[CH][C]=O(2676)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2407',
    isomers = [
        'C=C([O])[CH][C]=O(2676)',
    ],
    reactants = [
        ('HCCO(2227)', 'ketene(1375)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2407',
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


species(
    label = '[O][C]=CC=[C][O](9630)',
    structure = SMILES('[O][C]=C[CH][C]=O'),
    E0 = (401.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,533.66,533.804],'cm^-1')),
        HinderedRotor(inertia=(0.128472,'amu*angstrom^2'), symmetry=1, barrier=(25.9769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12849,'amu*angstrom^2'), symmetry=1, barrier=(25.9754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20645,0.0272905,8.57194e-06,-3.57664e-08,1.62277e-11,48348.4,23.0134], Tmin=(100,'K'), Tmax=(1001.37,'K')), NASAPolynomial(coeffs=[14.7008,0.00621477,-3.04802e-06,7.23519e-10,-6.08836e-14,44400.5,-44.4997], Tmin=(1001.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCJC=O) + radical(CCCJ=O) + radical(C=CJO)"""),
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
    label = '[O][C]=CC=C=O(9628)',
    structure = SMILES('O=[C]C=C[C]=O'),
    E0 = (116.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1850,1860,440,470,900,1000],'cm^-1')),
        HinderedRotor(inertia=(0.97355,'amu*angstrom^2'), symmetry=1, barrier=(22.3838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.971696,'amu*angstrom^2'), symmetry=1, barrier=(22.3412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28772,0.069745,-0.000137169,1.34309e-07,-5.00931e-11,14089.5,19.5988], Tmin=(100,'K'), Tmax=(809.64,'K')), NASAPolynomial(coeffs=[6.69836,0.0232224,-1.4311e-05,2.95401e-09,-2.10694e-13,13862,-1.3559], Tmin=(809.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CCJ=O)"""),
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
    label = '[O][C]=C[C][C]=O(10416)',
    structure = SMILES('[O][C]=[C][CH][C]=O'),
    E0 = (639.207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,562.219,562.433],'cm^-1')),
        HinderedRotor(inertia=(0.111242,'amu*angstrom^2'), symmetry=1, barrier=(24.9558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111173,'amu*angstrom^2'), symmetry=1, barrier=(24.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16909,0.0319165,-1.61792e-05,-6.7511e-09,5.46949e-12,76951.9,24.3159], Tmin=(100,'K'), Tmax=(1050.56,'K')), NASAPolynomial(coeffs=[13.6932,0.00541843,-3.15997e-06,7.34275e-10,-5.90938e-14,73571.5,-36.4143], Tmin=(1050.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCJC=O) + radical(Cds_S) + radical(CCCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[O]C1=CC1[C]=O(10417)',
    structure = SMILES('O=[C]C1[CH]C1=O'),
    E0 = (248.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75952,0.0147239,3.65219e-05,-6.26321e-08,2.63027e-11,29925.5,20.3976], Tmin=(100,'K'), Tmax=(939.774,'K')), NASAPolynomial(coeffs=[12.2009,0.00643157,-1.14847e-06,2.03071e-10,-1.93705e-14,26742.6,-32.0581], Tmin=(939.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[O]C=C=C[C]=O(10418)',
    structure = SMILES('[O]C=C=C[C]=O'),
    E0 = (146.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7611,0.0442043,-4.36984e-05,2.00881e-08,-3.57078e-12,17714.4,19.2867], Tmin=(100,'K'), Tmax=(1370.85,'K')), NASAPolynomial(coeffs=[14.1522,0.00804831,-4.13608e-06,8.48302e-10,-6.20339e-14,14317.1,-44.4016], Tmin=(1370.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=C=C[C]=O(10419)',
    structure = SMILES('O=[C][CH][C]=C=O'),
    E0 = (341.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,1855,455,950,2120,512.5,787.5,479.404],'cm^-1')),
        HinderedRotor(inertia=(0.215361,'amu*angstrom^2'), symmetry=1, barrier=(35.1958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215468,'amu*angstrom^2'), symmetry=1, barrier=(35.1993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22015,0.0335617,-2.53176e-05,3.89688e-09,2.09508e-12,41163.4,21.026], Tmin=(100,'K'), Tmax=(993.155,'K')), NASAPolynomial(coeffs=[11.7864,0.00675201,-2.5259e-06,4.78945e-10,-3.53408e-14,38685.3,-27.9698], Tmin=(993.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + missing(Cdd-CdO2d) + radical(CCJC=O) + radical(CCCJ=C=O) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][C]C[C]=O(9876)',
    structure = SMILES('O=[C][C]C[C]=O'),
    E0 = (458.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1850,1860,440,470,900,1000,182.144,182.677],'cm^-1')),
        HinderedRotor(inertia=(0.292823,'amu*angstrom^2'), symmetry=1, barrier=(6.91414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292429,'amu*angstrom^2'), symmetry=1, barrier=(6.91303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293143,'amu*angstrom^2'), symmetry=1, barrier=(6.91327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70761,0.0592809,-0.000112815,1.05872e-07,-3.68654e-11,55186.9,23.4387], Tmin=(100,'K'), Tmax=(888.668,'K')), NASAPolynomial(coeffs=[6.03776,0.0188705,-9.29494e-06,1.72269e-09,-1.13966e-13,55243.3,7.70708], Tmin=(888.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=[C][CH][C]=O(10420)',
    structure = SMILES('[O]C=[C][CH][C]=O'),
    E0 = (399.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,280.37,281.055],'cm^-1')),
        HinderedRotor(inertia=(0.536924,'amu*angstrom^2'), symmetry=1, barrier=(29.2866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519015,'amu*angstrom^2'), symmetry=1, barrier=(29.3345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96492,0.0301427,1.04884e-05,-4.48047e-08,2.11533e-11,48130.7,21.9171], Tmin=(100,'K'), Tmax=(973.854,'K')), NASAPolynomial(coeffs=[17.5085,0.00182392,-6.10857e-07,2.51427e-10,-2.89905e-14,43418.6,-61.309], Tmin=(973.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCJC=O) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=[C][CH]C=O(10421)',
    structure = SMILES('[O][C]=[C]C=C[O]'),
    E0 = (398.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.64757,'amu*angstrom^2'), symmetry=1, barrier=(37.8808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.551097,0.0545795,-6.13597e-05,3.11544e-08,-5.66199e-12,48011,25.6497], Tmin=(100,'K'), Tmax=(1652.03,'K')), NASAPolynomial(coeffs=[16.7817,-0.00067851,3.30396e-06,-7.87921e-10,5.6733e-14,44826.1,-54.2102], Tmin=(1652.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH][C]=[C]O(10422)',
    structure = SMILES('O=[C][CH][C]=[C]O'),
    E0 = (497.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,346.469],'cm^-1')),
        HinderedRotor(inertia=(0.260364,'amu*angstrom^2'), symmetry=1, barrier=(22.153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260317,'amu*angstrom^2'), symmetry=1, barrier=(22.1565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.260333,'amu*angstrom^2'), symmetry=1, barrier=(22.1556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75129,0.0390167,-2.03285e-05,-1.15227e-08,9.36186e-12,59955.1,24.402], Tmin=(100,'K'), Tmax=(980.563,'K')), NASAPolynomial(coeffs=[16.7011,0.00212716,-7.56006e-07,2.2986e-10,-2.35795e-14,55864.9,-53.3351], Tmin=(980.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Cds_S) + radical(C=CJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1=CC=C1[O](10423)',
    structure = SMILES('O=C1[CH][CH]C1=O'),
    E0 = (205.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84024,0.0182007,1.39686e-05,-2.75597e-08,1.04757e-11,24782.4,19.2218], Tmin=(100,'K'), Tmax=(1056.48,'K')), NASAPolynomial(coeffs=[8.24496,0.0153256,-6.92075e-06,1.37953e-09,-1.01096e-13,22658.9,-11.7951], Tmin=(1056.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclobutane) + radical(CCJCC=O) + radical(CCJCC=O)"""),
)

species(
    label = '[O][C]=[C]C=[C]O(10424)',
    structure = SMILES('[O][C]=[C]C=[C]O'),
    E0 = (496.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29509,'amu*angstrom^2'), symmetry=1, barrier=(29.7767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29046,'amu*angstrom^2'), symmetry=1, barrier=(29.6702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83242,0.0579291,-7.42006e-05,4.28558e-08,-8.9338e-12,59813.2,26.3366], Tmin=(100,'K'), Tmax=(1401.78,'K')), NASAPolynomial(coeffs=[17.0767,-0.00162875,3.66015e-06,-8.9363e-10,6.71606e-14,56556.4,-52.8918], Tmin=(1401.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH][C][O](10218)',
    structure = SMILES('[CH][C][O]'),
    E0 = (885.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([418.405,418.464,1386.35,1386.4,1890.01],'cm^-1')),
        HinderedRotor(inertia=(0.0607177,'amu*angstrom^2'), symmetry=1, barrier=(7.54278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.35351,0.0168371,-3.18729e-05,3.43969e-08,-1.43789e-11,106493,11.3039], Tmin=(100,'K'), Tmax=(740.605,'K')), NASAPolynomial(coeffs=[3.99585,0.00842558,-4.82643e-06,1.03992e-09,-7.72054e-14,106534,9.313], Tmin=(740.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(885.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCJ2_triplet) + radical(CH2_triplet)"""),
)

species(
    label = '[O][C][CH][C]=C[O](10425)',
    structure = SMILES('[O][CH][C][CH][C]=O'),
    E0 = (794.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,1855,455,950,180,180,180,1371.62],'cm^-1')),
        HinderedRotor(inertia=(0.164082,'amu*angstrom^2'), symmetry=1, barrier=(3.77258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163976,'amu*angstrom^2'), symmetry=1, barrier=(3.77012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0442551,'amu*angstrom^2'), symmetry=1, barrier=(59.1337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56468,0.0641788,-0.000125625,1.19849e-07,-4.20405e-11,95667.3,23.9752], Tmin=(100,'K'), Tmax=(894.767,'K')), NASAPolynomial(coeffs=[5.44966,0.0211271,-1.03955e-05,1.91333e-09,-1.25604e-13,96000.2,11.4097], Tmin=(894.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(794.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCHO) + radical(CCsJOH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]C1[CH]C1=O(10426)',
    structure = SMILES('[O][C]C1[CH]C1=O'),
    E0 = (661.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88727,0.0420824,-4.1278e-05,1.99494e-08,-3.77354e-12,79610,18.9124], Tmin=(100,'K'), Tmax=(1288.16,'K')), NASAPolynomial(coeffs=[12.1274,0.0102848,-4.25152e-06,7.86996e-10,-5.46252e-14,76971.7,-33.0834], Tmin=(1288.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CCOJ) + radical(CCJC=O) + radical(CH2_triplet)"""),
)

species(
    label = '[O][C]C1C=[C]O1(10427)',
    structure = SMILES('[O][C]C1C=[C]O1'),
    E0 = (729.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94669,0.0335553,-4.50039e-06,-2.77319e-08,1.52421e-11,87776.5,20.5722], Tmin=(100,'K'), Tmax=(959.672,'K')), NASAPolynomial(coeffs=[16.1492,0.00264288,-3.932e-07,1.26829e-10,-1.58712e-14,83748.1,-54.1479], Tmin=(959.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CCOJ) + radical(C=CJO) + radical(CH2_triplet)"""),
)

species(
    label = 'O=C1[CH][CH][C]O1(10428)',
    structure = SMILES('[O]C1=C[CH][C]O1'),
    E0 = (433.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94748,0.0305601,1.51774e-05,-5.73917e-08,2.91036e-11,52250.6,10.8366], Tmin=(100,'K'), Tmax=(901.74,'K')), NASAPolynomial(coeffs=[18.1211,-0.00195119,3.99682e-06,-8.77236e-10,5.88848e-14,47738.6,-74.3632], Tmin=(901.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + ring(2,3-Dihydrofuran) + radical(C=COJ) + radical(C=CCJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[C]1[CH]C=[C]OO1(10429)',
    structure = SMILES('[C]1[CH]C=[C]OO1'),
    E0 = (779.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13322,0.0319958,-3.30088e-06,-2.59129e-08,1.44888e-11,93872.1,14.1086], Tmin=(100,'K'), Tmax=(923.574,'K')), NASAPolynomial(coeffs=[13.2325,0.00736729,-1.37498e-06,1.79892e-10,-1.34522e-14,90822.1,-43.9695], Tmin=(923.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(779.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(C=CCJCO) + radical(CH2_triplet) + radical(C=CJO)"""),
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
    label = '[CH]=C[C][O](10430)',
    structure = SMILES('[CH][CH][C]=O'),
    E0 = (571.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,180,1371.68,1374.81],'cm^-1')),
        HinderedRotor(inertia=(8.95207e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383504,'amu*angstrom^2'), symmetry=1, barrier=(8.81751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8972,0.0246759,-2.89293e-05,1.86516e-08,-4.76072e-12,68758.2,15.1597], Tmin=(100,'K'), Tmax=(960.205,'K')), NASAPolynomial(coeffs=[6.91501,0.00793846,-2.78248e-06,4.97856e-10,-3.41507e-14,67986.6,-4.06073], Tmin=(960.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJCHO) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = 'O=[C][CH]C1[C]O1(10431)',
    structure = SMILES('O=[C][CH]C1[C]O1'),
    E0 = (610.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50015,0.0491059,-6.35712e-05,4.00176e-08,-9.15205e-12,73562.4,22.0398], Tmin=(100,'K'), Tmax=(1299.96,'K')), NASAPolynomial(coeffs=[11.8467,0.00501125,1.45325e-06,-5.83061e-10,5.10254e-14,71908.1,-26.6069], Tmin=(1299.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCJCO) + radical(CH2_triplet) + radical(CCCJ=O)"""),
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
    E0 = (401.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (401.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (905.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (850.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (404.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (424.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (561.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (607.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (640.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (580.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (548.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (625.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (409.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (644.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1105.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (817.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (663.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (730.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (475.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (780.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1065.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (613.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (472.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['HCCO(2227)', 'HCCO(2227)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[O][C]=CC=C=O(9628)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C][O](6861)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[O][C]=C[C][C]=O(10416)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[O]C1=CC1[C]=O(10417)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.36e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[O]C=C=C[C]=O(10418)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[O][C]=C=C[C]=O(10419)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCCO(2227)', '[CH]=[C][O](6861)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(72.9547,'m^3/(mol*s)'), n=1.66457, Ea=(16.7701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;YJ] for rate rule [Ct_Ct;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=[C][C]C[C]=O(9876)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.28974e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=[C][CH][C]=O(10420)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.97782e+08,'s^-1'), n=1.48417, Ea=(181.505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=[C][CH]C=O(10421)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(781966,'s^-1'), n=1.9774, Ea=(150.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out;XH_out] for rate rule [R3HJ;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=[C][CH][C]=[C]O(10422)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2113.06,'s^-1'), n=2.4992, Ea=(127.803,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3Hall;Y_rad_out;O_H_out] + [R3Hall;Cd_rad_out;XH_out] for rate rule [R3HJ;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[O]C1=CC=C1[O](10423)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=[C]C=[C]O(10424)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1194.19,'s^-1'), n=2.56519, Ea=(148.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;O_H_out] + [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HCCO(2227)', '[CH][C][O](10218)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C][CH][C]=C[O](10425)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[O][C]C1[CH]C1=O(10426)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.22706e+11,'s^-1'), n=0.395207, Ea=(261.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[O][C]C1C=[C]O1(10427)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.74701e+10,'s^-1'), n=0.685238, Ea=(328.663,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['O=C1[CH][CH][C]O1(10428)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.10072e+11,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['[C]1[CH]C=[C]OO1(10429)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.41956e+10,'s^-1'), n=0.267163, Ea=(379.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;multiplebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]=O(2355)', '[CH]=C[C][O](10430)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][C]=CC=[C][O](9630)'],
    products = ['O=[C][CH]C1[C]O1(10431)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.61926e+11,'s^-1'), n=0.419784, Ea=(211.797,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R4;multiplebond_intra;radadd_intra_O] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CO(2039)', '[CH]=C[C][O](10430)'],
    products = ['[O][C]=CC=[C][O](9630)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

network(
    label = 'PDepNetwork #2552',
    isomers = [
        '[O][C]=CC=[C][O](9630)',
    ],
    reactants = [
        ('HCCO(2227)', 'HCCO(2227)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2552',
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


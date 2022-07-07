species(
    label = '[CH]=CC([O])=C[O](15072)',
    structure = SMILES('[CH]=CC([O])=C[O]'),
    E0 = (197.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2425,'amu*angstrom^2'), symmetry=1, barrier=(28.5675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270776,0.0600526,-6.76235e-05,3.43253e-08,-6.25885e-12,23918.9,23.8081], Tmin=(100,'K'), Tmax=(1632.35,'K')), NASAPolynomial(coeffs=[18.5321,-0.000982274,3.42858e-06,-8.05211e-10,5.7624e-14,20126.9,-66.5945], Tmin=(1632.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08745,0.00578572,8.56332e-06,-1.72824e-08,7.83594e-12,26273.8,4.4608], Tmin=(100,'K'), Tmax=(897.945,'K')), NASAPolynomial(coeffs=[5.89068,0.00209055,4.88462e-08,-5.66903e-11,4.15085e-15,25415.9,-10.7352], Tmin=(897.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.303,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]C=[C]C=O(15283)',
    structure = SMILES('[CH]C=C=C[O]'),
    E0 = (450.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07021,'amu*angstrom^2'), symmetry=1, barrier=(47.5981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08056,0.032167,3.59315e-06,-3.27881e-08,1.63885e-11,54207.6,17.2931], Tmin=(100,'K'), Tmax=(940.198,'K')), NASAPolynomial(coeffs=[12.873,0.0113851,-3.34975e-06,5.67578e-10,-4.11607e-14,51067.4,-40.0163], Tmin=(940.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=[CH])[O](15067)',
    structure = SMILES('[CH]=CC(=[CH])[O]'),
    E0 = (512.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3115,3125,620,680,785,800,1600,1700],'cm^-1')),
        HinderedRotor(inertia=(1.01991,'amu*angstrom^2'), symmetry=1, barrier=(23.4498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3512.15,'J/mol'), sigma=(5.74623,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=548.59 K, Pc=42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55548,0.0463587,-5.45809e-05,3.0561e-08,-6.36527e-12,61676.8,18.3157], Tmin=(100,'K'), Tmax=(1338.62,'K')), NASAPolynomial(coeffs=[13.8308,0.00354016,2.7773e-07,-1.85393e-10,1.66465e-14,58940.4,-42.4316], Tmin=(1338.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[C]=CC([O])=C[O](15422)',
    structure = SMILES('[C]=CC([O])=C[O]'),
    E0 = (508.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18236,'amu*angstrom^2'), symmetry=1, barrier=(27.1847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770767,0.0576609,-7.06945e-05,3.87987e-08,-7.72628e-12,61298.9,21.9734], Tmin=(100,'K'), Tmax=(1448.42,'K')), NASAPolynomial(coeffs=[18.1565,-0.0019868,3.12636e-06,-7.2484e-10,5.27734e-14,57482.9,-64.1306], Tmin=(1448.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[O]C=C1C=CO1(15416)',
    structure = SMILES('[O]C=C1C=CO1'),
    E0 = (29.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93703,0.0229519,5.46217e-05,-1.06602e-07,4.8071e-11,3694.11,15.9933], Tmin=(100,'K'), Tmax=(914.186,'K')), NASAPolynomial(coeffs=[22.0244,-0.00630158,6.40706e-06,-1.27802e-09,8.08365e-14,-2428.9,-92.5158], Tmin=(914.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC1=COO1(15423)',
    structure = SMILES('[CH]=CC1=COO1'),
    E0 = (437.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92864,0.0331057,4.26415e-06,-3.76732e-08,1.88921e-11,52726.8,18.1452], Tmin=(100,'K'), Tmax=(951.449,'K')), NASAPolynomial(coeffs=[15.7702,0.00575406,-1.23492e-06,2.474e-10,-2.3234e-14,48697,-55.2789], Tmin=(951.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P)"""),
)

species(
    label = '[O]C1=COC=C1(15424)',
    structure = SMILES('[O]C1=COC=C1'),
    E0 = (-123.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54514,0.0172316,4.1722e-05,-7.21951e-08,3.0343e-11,-14762.1,15.9579], Tmin=(100,'K'), Tmax=(942.837,'K')), NASAPolynomial(coeffs=[13.8401,0.00674189,-1.13759e-06,2.1592e-10,-2.18844e-14,-18555.5,-46.6915], Tmin=(942.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Furan) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC(O)=C=O(15425)',
    structure = SMILES('[CH]=CC(O)=C=O'),
    E0 = (114.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27484,0.0568111,-6.81771e-05,3.6252e-08,-6.19417e-12,13921,19.2855], Tmin=(100,'K'), Tmax=(888.032,'K')), NASAPolynomial(coeffs=[15.0304,0.0060076,-1.20774e-06,1.23377e-10,-5.81834e-15,11038.1,-47.9204], Tmin=(888.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsOs) + group(Cd-Cd(CCO)H) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(Cds_P)"""),
)

species(
    label = 'C=CC([O])=C=O(15426)',
    structure = SMILES('C=CC(=O)[C]=O'),
    E0 = (-22.5993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32584,0.0399421,-4.03286e-05,2.34923e-08,-5.88854e-12,-2660.41,19.213], Tmin=(100,'K'), Tmax=(935.076,'K')), NASAPolynomial(coeffs=[6.7188,0.0211503,-1.01839e-05,2.00063e-09,-1.42598e-13,-3481.96,-1.68563], Tmin=(935.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.5993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C]C([O])[CH][O](15427)',
    structure = SMILES('[CH]=[C]C([O])[CH][O]'),
    E0 = (772.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,343.199,343.3,1354.61],'cm^-1')),
        HinderedRotor(inertia=(0.0878627,'amu*angstrom^2'), symmetry=1, barrier=(7.30651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868527,'amu*angstrom^2'), symmetry=1, barrier=(7.31616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46724,0.0636721,-0.00011299,1.04362e-07,-3.67754e-11,92951.2,26.2332], Tmin=(100,'K'), Tmax=(853.62,'K')), NASAPolynomial(coeffs=[6.66929,0.0226093,-1.15119e-05,2.20928e-09,-1.51114e-13,92671.1,5.52067], Tmin=(853.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C]1OC1[O](15428)',
    structure = SMILES('[CH]C=C1OC1[O]'),
    E0 = (377.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83078,0.039085,-1.20239e-05,-1.39253e-08,8.36878e-12,45515.7,18.4494], Tmin=(100,'K'), Tmax=(1003.82,'K')), NASAPolynomial(coeffs=[12.271,0.0164371,-6.50411e-06,1.21859e-09,-8.72814e-14,42464.7,-36.715], Tmin=(1003.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(methyleneoxirane) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC1([O])[CH]O1(15429)',
    structure = SMILES('[CH]=CC1([O])[CH]O1'),
    E0 = (451.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.13202,0.0646685,-7.87797e-05,4.36211e-08,-8.43879e-12,54430.4,24.9423], Tmin=(100,'K'), Tmax=(1610.07,'K')), NASAPolynomial(coeffs=[16.2522,-0.00159516,6.76587e-06,-1.65955e-09,1.22988e-13,52467.4,-51.6173], Tmin=(1610.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CC(C)(O)OJ) + radical(CCsJO) + radical(Cds_P)"""),
)

species(
    label = '[O][C]1C=CC1[O](15430)',
    structure = SMILES('[O]C1=C[CH]C1[O]'),
    E0 = (242.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00495,0.0259021,3.5937e-05,-7.7434e-08,3.4742e-11,29216.5,17.6007], Tmin=(100,'K'), Tmax=(934.044,'K')), NASAPolynomial(coeffs=[18.4499,0.00138474,1.58643e-06,-2.9703e-10,1.22472e-14,24141.9,-71.3343], Tmin=(934.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O][CH]C1([O])C=C1(15431)',
    structure = SMILES('[O][CH]C1([O])C=C1'),
    E0 = (520.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15378,0.0447948,-4.98527e-05,3.19205e-08,-8.75278e-12,62681.4,18.881], Tmin=(100,'K'), Tmax=(861.753,'K')), NASAPolynomial(coeffs=[6.88764,0.0228209,-1.16029e-05,2.32887e-09,-1.67796e-13,61865.6,-3.25264], Tmin=(861.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C#CC([O])=C[O](13171)',
    structure = SMILES('C#CC([O])=C[O]'),
    E0 = (127.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.58154,'amu*angstrom^2'), symmetry=1, barrier=(36.3626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605988,0.0547459,-6.14813e-05,3.08789e-08,-5.58321e-12,15448.7,21.0975], Tmin=(100,'K'), Tmax=(1636.22,'K')), NASAPolynomial(coeffs=[17.7517,-0.00135695,2.95693e-06,-6.75319e-10,4.76604e-14,11737,-64.2598], Tmin=(1636.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC([O])=C=O(15432)',
    structure = SMILES('[CH]=CC(=O)[C]=O'),
    E0 = (224.497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,3120,650,792.5,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.522709,'amu*angstrom^2'), symmetry=1, barrier=(12.0181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1168,'amu*angstrom^2'), symmetry=1, barrier=(25.6775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17008,0.0443164,-5.91984e-05,4.20504e-08,-1.11324e-11,27062.8,20.1808], Tmin=(100,'K'), Tmax=(638.098,'K')), NASAPolynomial(coeffs=[7.16358,0.0179795,-8.95992e-06,1.75792e-09,-1.2415e-13,26324.5,-2.45865], Tmin=(638.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_P) + radical(CCCJ=O)"""),
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
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]C([O])=C[O](13159)',
    structure = SMILES('[CH]=[C]C([O])=C[O]'),
    E0 = (396.595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55203,'amu*angstrom^2'), symmetry=1, barrier=(35.6842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690677,0.0610483,-7.87226e-05,4.57036e-08,-9.59389e-12,47828.7,21.8755], Tmin=(100,'K'), Tmax=(1385.71,'K')), NASAPolynomial(coeffs=[17.8401,-0.00178865,3.73025e-06,-9.08451e-10,6.82983e-14,44356,-61.8356], Tmin=(1385.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([O])=[C][O](15433)',
    structure = SMILES('[CH]=CC([O])=[C][O]'),
    E0 = (437.343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3120,650,792.5,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.979774,'amu*angstrom^2'), symmetry=1, barrier=(22.5269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19025,0.0536825,-6.70689e-05,3.87117e-08,-8.2425e-12,52708.5,23.6199], Tmin=(100,'K'), Tmax=(1315.64,'K')), NASAPolynomial(coeffs=[15.8793,0.00134649,1.35282e-06,-3.94088e-10,3.11572e-14,49507.7,-48.7512], Tmin=(1315.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C([O])=C[O](11900)',
    structure = SMILES('C=[C]C([O])=C[O]'),
    E0 = (149.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55298,'amu*angstrom^2'), symmetry=1, barrier=(35.706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562656,0.0598018,-6.95646e-05,3.74263e-08,-7.30439e-12,18117.9,21.9387], Tmin=(100,'K'), Tmax=(1502.94,'K')), NASAPolynomial(coeffs=[17.396,0.00121988,2.6563e-06,-7.09798e-10,5.39892e-14,14614.4,-60.9531], Tmin=(1502.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC([O])=[C]O(15434)',
    structure = SMILES('[CH]=CC([O])=[C]O'),
    E0 = (295.881,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.923153,'amu*angstrom^2'), symmetry=1, barrier=(21.2251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928946,'amu*angstrom^2'), symmetry=1, barrier=(21.3583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544173,0.0634832,-8.06967e-05,4.62669e-08,-9.6109e-12,35721.5,24.5244], Tmin=(100,'K'), Tmax=(1394.16,'K')), NASAPolynomial(coeffs=[18.6433,-0.00168628,3.66644e-06,-8.86636e-10,6.6248e-14,31961.7,-64.1929], Tmin=(1394.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C(O)=C[O](15435)',
    structure = SMILES('[CH]=[C]C(O)=C[O]'),
    E0 = (258.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.41888,'amu*angstrom^2'), symmetry=1, barrier=(32.6228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41922,'amu*angstrom^2'), symmetry=1, barrier=(32.6307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0439361,0.0692845,-8.57578e-05,4.70221e-08,-9.22648e-12,31287.8,23.0163], Tmin=(100,'K'), Tmax=(1509.54,'K')), NASAPolynomial(coeffs=[20.7967,-0.00417424,5.3564e-06,-1.21943e-09,8.82456e-14,27073.5,-79.2284], Tmin=(1509.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(O)=[C][O](15436)',
    structure = SMILES('[CH]C=C(O)[C]=O'),
    E0 = (207.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.94723,'amu*angstrom^2'), symmetry=1, barrier=(44.7707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94673,'amu*angstrom^2'), symmetry=1, barrier=(44.7592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.94725,'amu*angstrom^2'), symmetry=1, barrier=(44.7712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954089,0.0658013,-7.84806e-05,4.61092e-08,-1.05785e-11,25053.4,20.4795], Tmin=(100,'K'), Tmax=(1069.28,'K')), NASAPolynomial(coeffs=[14.6695,0.014494,-6.50582e-06,1.23462e-09,-8.66755e-14,22120.3,-46.608], Tmin=(1069.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC([O])=[C][O](15437)',
    structure = SMILES('C=CC([O])=[C][O]'),
    E0 = (190.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07559,0.0522916,-5.74638e-05,2.99306e-08,-5.76913e-12,22997.1,23.6343], Tmin=(100,'K'), Tmax=(1473.45,'K')), NASAPolynomial(coeffs=[15.7197,0.00395646,4.77503e-07,-2.37412e-10,2.00408e-14,19613.1,-49.5303], Tmin=(1473.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]C([O])=CO(15438)',
    structure = SMILES('[CH]=C=C([O])[CH]O'),
    E0 = (188.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979579,0.0516709,-5.57116e-05,2.84921e-08,-5.34039e-12,22802.9,25.4059], Tmin=(100,'K'), Tmax=(1554.8,'K')), NASAPolynomial(coeffs=[15.3126,0.00383137,1.02073e-06,-3.69647e-10,2.97185e-14,19671.3,-45.806], Tmin=(1554.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=C=CJ)"""),
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
    label = '[CH]=C[CH][O](15182)',
    structure = SMILES('[CH]C=C[O]'),
    E0 = (345.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,278.044,280.154,284.215,284.973],'cm^-1')),
        HinderedRotor(inertia=(0.978465,'amu*angstrom^2'), symmetry=1, barrier=(51.1711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06273,0.0225089,-7.83191e-06,1.77811e-11,2.905e-13,41549.3,13.3388], Tmin=(100,'K'), Tmax=(1878.6,'K')), NASAPolynomial(coeffs=[9.20694,0.015113,-6.46712e-06,1.14479e-09,-7.39143e-14,38237.4,-22.8483], Tmin=(1878.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([O])[C]=O(15439)',
    structure = SMILES('[CH][CH]C([O])[C]=O'),
    E0 = (654.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1855,455,950,180,790.563,790.605,793.543,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00721914,'amu*angstrom^2'), symmetry=1, barrier=(3.21123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0071864,'amu*angstrom^2'), symmetry=1, barrier=(3.19403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236733,'amu*angstrom^2'), symmetry=1, barrier=(10.4305,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89298,0.0410922,-3.48601e-05,1.07444e-08,2.43474e-13,78837.9,29.1113], Tmin=(100,'K'), Tmax=(981.676,'K')), NASAPolynomial(coeffs=[12.3024,0.0095495,-3.27581e-06,5.77139e-10,-4.04117e-14,76270.3,-23.5841], Tmin=(981.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]=C([O])C[O](15440)',
    structure = SMILES('[CH][C]=C([O])C[O]'),
    E0 = (574.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,497.11,497.115,497.12,497.122,497.126,497.126],'cm^-1')),
        HinderedRotor(inertia=(0.00068213,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.304376,'amu*angstrom^2'), symmetry=1, barrier=(53.3806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83212,0.0546755,-8.35167e-05,7.75718e-08,-2.81119e-11,69160.5,24.5091], Tmin=(100,'K'), Tmax=(861.435,'K')), NASAPolynomial(coeffs=[3.5832,0.030721,-1.42519e-05,2.64402e-09,-1.7843e-13,69445.9,19.7301], Tmin=(861.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=C([O])[CH]O(15441)',
    structure = SMILES('[CH][C]=C([O])[CH]O'),
    E0 = (466.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,1685,370,521.57,521.571,521.571,521.571,521.573],'cm^-1')),
        HinderedRotor(inertia=(0.268696,'amu*angstrom^2'), symmetry=1, barrier=(51.8703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268699,'amu*angstrom^2'), symmetry=1, barrier=(51.8704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268697,'amu*angstrom^2'), symmetry=1, barrier=(51.8703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56327,0.0477141,-3.58927e-05,7.73849e-09,1.95755e-12,56143.5,25.2104], Tmin=(100,'K'), Tmax=(956.649,'K')), NASAPolynomial(coeffs=[12.6445,0.0147291,-5.10323e-06,8.67766e-10,-5.86172e-14,53412.5,-30.9516], Tmin=(956.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1O[C]1C=O(15442)',
    structure = SMILES('[CH]C1OC1=C[O]'),
    E0 = (362.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50755,0.0297532,4.77135e-05,-1.09264e-07,5.14244e-11,43754.4,18.4252], Tmin=(100,'K'), Tmax=(911.034,'K')), NASAPolynomial(coeffs=[26.3512,-0.0132651,9.77444e-06,-1.90837e-09,1.23198e-13,36486.3,-114.162], Tmin=(911.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C=C1[CH]OO1(15443)',
    structure = SMILES('[CH]C=C1[CH]OO1'),
    E0 = (542.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59458,0.0120022,7.41054e-05,-1.04174e-07,3.99264e-11,65353.7,20.733], Tmin=(100,'K'), Tmax=(973.125,'K')), NASAPolynomial(coeffs=[13.0812,0.0165613,-6.39277e-06,1.30666e-09,-1.02741e-13,61055.9,-41.169], Tmin=(973.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'HCO(1372)',
    structure = SMILES('[CH]=O'),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C[C]=O(15128)',
    structure = SMILES('[CH]C=C=O'),
    E0 = (292.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2120,512.5,787.5,246.578,247.126,247.139],'cm^-1')),
        HinderedRotor(inertia=(1.16013,'amu*angstrom^2'), symmetry=1, barrier=(50.7263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00497,0.0183622,-5.86157e-07,-8.19965e-09,3.38646e-12,35273.9,13.1075], Tmin=(100,'K'), Tmax=(1106.37,'K')), NASAPolynomial(coeffs=[5.96287,0.0150316,-6.05412e-06,1.11098e-09,-7.677e-14,34168.7,-3.49831], Tmin=(1106.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([C]=O)C=O(15444)',
    structure = SMILES('[CH]C([C]=O)C=O'),
    E0 = (291.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,606.205,606.205,606.205],'cm^-1')),
        HinderedRotor(inertia=(0.011995,'amu*angstrom^2'), symmetry=1, barrier=(3.12799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.011995,'amu*angstrom^2'), symmetry=1, barrier=(3.12799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.011995,'amu*angstrom^2'), symmetry=1, barrier=(3.12799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14414,0.0352903,-2.32583e-05,5.21285e-09,3.12755e-13,35135.6,23.9718], Tmin=(100,'K'), Tmax=(1152.52,'K')), NASAPolynomial(coeffs=[10.4612,0.0142421,-6.03846e-06,1.1374e-09,-7.98137e-14,32699.3,-19.5864], Tmin=(1152.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH][CH](6981)',
    structure = SMILES('[CH][CH]'),
    E0 = (800.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([355.88,360.972,1764.15,1768.16,1768.69],'cm^-1')),
        HinderedRotor(inertia=(0.000315567,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62974,0.00646269,3.64725e-07,-4.23355e-09,1.88225e-12,96340.6,7.14898], Tmin=(100,'K'), Tmax=(1051.09,'K')), NASAPolynomial(coeffs=[5.42018,0.00355693,-1.06544e-06,2.10792e-10,-1.61554e-14,95748.4,-2.60499], Tmin=(1051.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1OC=C1[O](13294)',
    structure = SMILES('[CH]C1OC=C1[O]'),
    E0 = (333.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62791,0.0266315,5.57588e-05,-1.16909e-07,5.40508e-11,40276,18.8551], Tmin=(100,'K'), Tmax=(909.925,'K')), NASAPolynomial(coeffs=[25.8995,-0.0127639,9.75515e-06,-1.91762e-09,1.24091e-13,33072.7,-111.26], Tmin=(909.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1C(=O)C1[O](15445)',
    structure = SMILES('[CH]C1C(=O)C1[O]'),
    E0 = (494.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09744,0.0295039,1.19686e-05,-4.577e-08,2.23308e-11,59533.3,21.3788], Tmin=(100,'K'), Tmax=(926.475,'K')), NASAPolynomial(coeffs=[15.086,0.00520539,-1.42167e-07,-3.27211e-11,-1.10041e-15,55762.7,-47.652], Tmin=(926.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(C=OCOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]CC(=O)[C]=O(15446)',
    structure = SMILES('[CH]CC(=O)[C]=O'),
    E0 = (296.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1855,455,950,389.351,389.654,389.84,2163.19],'cm^-1')),
        HinderedRotor(inertia=(0.00301103,'amu*angstrom^2'), symmetry=1, barrier=(9.99493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0926837,'amu*angstrom^2'), symmetry=1, barrier=(9.99112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.393443,'amu*angstrom^2'), symmetry=1, barrier=(42.3263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66856,0.0563363,-8.86804e-05,7.70178e-08,-2.64898e-11,35783.2,22.0156], Tmin=(100,'K'), Tmax=(820.106,'K')), NASAPolynomial(coeffs=[7.27986,0.0212753,-1.04829e-05,2.0135e-09,-1.39098e-13,35121.5,-2.36571], Tmin=(820.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC(=O)C=O(15071)',
    structure = SMILES('[CH]=CC(=O)C=O'),
    E0 = (64.5363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3557,0.0400806,-4.02297e-05,2.40833e-08,-6.38877e-12,7817.79,18.6942], Tmin=(100,'K'), Tmax=(873.879,'K')), NASAPolynomial(coeffs=[5.98345,0.0234753,-1.17267e-05,2.33877e-09,-1.68041e-13,7183.75,1.6815], Tmin=(873.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.5363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([O])[C]=O(15070)',
    structure = SMILES('[CH]=CC([O])[C]=O'),
    E0 = (380.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.159142,'amu*angstrom^2'), symmetry=1, barrier=(3.65898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15795,'amu*angstrom^2'), symmetry=1, barrier=(3.63157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98099,0.0419684,-4.26841e-05,2.22481e-08,-4.58545e-12,45832.6,25.5053], Tmin=(100,'K'), Tmax=(1180.53,'K')), NASAPolynomial(coeffs=[10.769,0.0121916,-4.84871e-06,8.81418e-10,-6.05747e-14,43757.8,-18.3499], Tmin=(1180.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    label = '[O]C1C=CC1=O(15447)',
    structure = SMILES('[O]C1C=CC1=O'),
    E0 = (105.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76167,0.0182301,2.14531e-05,-3.71783e-08,1.40091e-11,12750.2,18.3991], Tmin=(100,'K'), Tmax=(1043.52,'K')), NASAPolynomial(coeffs=[9.05533,0.0161805,-7.33259e-06,1.48438e-09,-1.10269e-13,10234.8,-17.9911], Tmin=(1043.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH][C]1C=CO1(15448)',
    structure = SMILES('[O]C=C1[CH][CH]O1'),
    E0 = (188.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24921,0.0357103,4.04372e-05,-1.09621e-07,5.45472e-11,22767,14.1223], Tmin=(100,'K'), Tmax=(888.335,'K')), NASAPolynomial(coeffs=[27.8245,-0.0159255,1.27589e-05,-2.64425e-09,1.81128e-13,15361.3,-126.049], Tmin=(888.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(C=CCJCO) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH]=[C]C(=O)C[O](15449)',
    structure = SMILES('[CH]=C=C([O])C[O]'),
    E0 = (296.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,540,610,2055,3120,650,792.5,1650,180,2023.55],'cm^-1')),
        HinderedRotor(inertia=(0.450879,'amu*angstrom^2'), symmetry=1, barrier=(10.3666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85663,0.0514917,-7.83838e-05,6.57365e-08,-2.1294e-11,35793.3,22.5193], Tmin=(100,'K'), Tmax=(906.294,'K')), NASAPolynomial(coeffs=[6.85135,0.0191728,-7.88808e-06,1.3713e-09,-8.84078e-14,35309.9,1.24202], Tmin=(906.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = '[O]C1(C=O)C=C1(15450)',
    structure = SMILES('[O]C1(C=O)C=C1'),
    E0 = (217.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25084,0.0412593,-4.76893e-05,3.3036e-08,-9.68286e-12,26262,18.6809], Tmin=(100,'K'), Tmax=(817.068,'K')), NASAPolynomial(coeffs=[6.50979,0.0204091,-9.41103e-06,1.80325e-09,-1.26342e-13,25566,-1.00555], Tmin=(817.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C=C([O])C=O(11930)',
    structure = SMILES('C=[C]C(=O)C=O'),
    E0 = (61.4013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,413.184],'cm^-1')),
        HinderedRotor(inertia=(0.569432,'amu*angstrom^2'), symmetry=1, barrier=(13.0924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03685,'amu*angstrom^2'), symmetry=1, barrier=(23.8392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34268,0.0408215,-4.11136e-05,2.37985e-08,-6.04675e-12,7440.87,19.0438], Tmin=(100,'K'), Tmax=(910.792,'K')), NASAPolynomial(coeffs=[6.38375,0.0230738,-1.18843e-05,2.40342e-09,-1.74035e-13,6704.76,-0.074362], Tmin=(910.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.4013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(C=CJC=O)"""),
)

species(
    label = '[O][C]1[CH]OC=C1(15451)',
    structure = SMILES('[O]C1=C[CH]O[CH]1'),
    E0 = (65.5394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82948,0.0160763,3.13305e-05,-4.86336e-08,1.86087e-11,7933.28,17.822], Tmin=(100,'K'), Tmax=(992.306,'K')), NASAPolynomial(coeffs=[8.43901,0.0173928,-6.83082e-06,1.30568e-09,-9.53173e-14,5641.91,-15.1335], Tmin=(992.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.5394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(C=CCJ(O)C) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=[C]C([O])C=O(15452)',
    structure = SMILES('[CH]=[C]C([O])C=O'),
    E0 = (458.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0873983,'amu*angstrom^2'), symmetry=1, barrier=(2.00946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0877962,'amu*angstrom^2'), symmetry=1, barrier=(2.01861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22089,0.0411605,-4.35438e-05,2.53377e-08,-6.10586e-12,55187.2,24.3058], Tmin=(100,'K'), Tmax=(990.052,'K')), NASAPolynomial(coeffs=[7.97788,0.017901,-8.30381e-06,1.6082e-09,-1.13832e-13,54047.2,-3.41077], Tmin=(990.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    E0 = (197.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (969.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1031.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (720.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (205.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (437.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (204.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (260.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (222.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (795.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (384.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (451.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (335.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (520.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (358.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (474.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (561.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (450.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (795.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (608.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (649.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (365.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (458.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (409.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (400.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (601.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (311.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (531.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (677.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (637.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (491.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (384.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (542.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (419.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (603.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (414.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (780.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (439.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (494.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (426.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (197.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (543.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (779.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (205.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (335.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (438.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (217.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (225.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (238.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (640.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['OCHCO(3676)', 'C2H2(1342)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH]C=[C]C=O(15283)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]=CC(=[CH])[O](15067)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[C]=CC([O])=C[O](15422)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O]C=C1C=CO1(15416)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]=CC1=COO1(15423)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O]C1=COC=C1(15424)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]=CC(O)=C=O(15425)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['C=CC([O])=C=O(15426)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]C([O])[CH][O](15427)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]=C[C]1OC1[O](15428)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]=CC1([O])[CH]O1(15429)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(253.718,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O][C]1C=CC1[O](15430)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O][CH]C1([O])C=C1(15431)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(323.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 322.6 to 323.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C#CC([O])=C[O](13171)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.642e+08,'cm^3/(mol*s)'), n=1.548, Ea=(19.0205,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 195 used for Ct-Cd_Ct-H;HJ
Exact match found for rate rule [Ct-Cd_Ct-H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]=CC([O])=C=O(15432)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OCHCO(3676)', 'C2H2(T)(1343)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7985.59,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C]=C[O](9592)', 'C2H2(1342)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(46.4627,'m^3/(mol*s)'), n=1.51997, Ea=(27.4714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-H;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=C[O](9592)', 'C2H2(T)(1343)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.62e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH]=[C]C([O])=C[O](13159)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH]=CC([O])=[C][O](15433)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['C=[C]C([O])=C[O](11900)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=CC([O])=[C]O(15434)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]C(O)=C[O](15435)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.87605e+06,'s^-1'), n=1.79171, Ea=(150.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=CC(O)=[C][O](15436)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['C=CC([O])=[C][O](15437)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]=[C]C([O])=CO(15438)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CO(2039)', '[CH]=C[CH][O](15182)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH][CH]C([O])[C]=O(15439)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH][C]=C([O])C[O](15440)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH][C]=C([O])[CH]O(15441)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]C1O[C]1C=O(15442)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]C=C1[CH]OO1(15443)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(345.22,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 340.8 to 345.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['HCO(1372)', '[CH]=C[C]=O(15128)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;CO_pri_rad] for rate rule [Ck_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['HCO(1372)', '[CH]=C[C][O](10430)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]C([C]=O)C=O(15444)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['OCHCO(3676)', '[CH][CH](6981)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]C1OC=C1[O](13294)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]C1C(=O)C1[O](15445)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.63856e+09,'s^-1'), n=0.755479, Ea=(296.724,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 294.0 to 296.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]CC(=O)[C]=O(15446)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(11891.5,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;XH_out] for rate rule [R3H_SS;CO_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[CH]=CC(=O)C=O(15071)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=CC([O])[C]=O(15070)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH][O](1548)', '[CH]=C[C]=O(15128)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O]C1C=CC1=O(15447)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O][CH][C]1C=CO1(15448)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=[C]C(=O)C[O](15449)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O]C1(C=O)C=C1(15450)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(20.2523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination
Ea raised from 18.9 to 20.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['C=C=C([O])C=O(11930)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=CC([O])=C[O](15072)'],
    products = ['[O][C]1[CH]OC=C1(15451)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.95882e+10,'s^-1'), n=0.283532, Ea=(41.0907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=[C]C([O])C=O(15452)'],
    products = ['[CH]=CC([O])=C[O](15072)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(3.64487e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3484',
    isomers = [
        '[CH]=CC([O])=C[O](15072)',
    ],
    reactants = [
        ('OCHCO(3676)', 'C2H2(1342)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3484',
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


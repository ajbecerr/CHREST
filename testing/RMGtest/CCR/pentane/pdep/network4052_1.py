species(
    label = 'C=[C][CH][CH][C]=O(16165)',
    structure = SMILES('[CH2][C]=C[CH][C]=O'),
    E0 = (582.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,1855,455,950,612.29],'cm^-1')),
        HinderedRotor(inertia=(0.095959,'amu*angstrom^2'), symmetry=1, barrier=(25.5686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0959819,'amu*angstrom^2'), symmetry=1, barrier=(25.5409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961329,'amu*angstrom^2'), symmetry=1, barrier=(25.5573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92363,0.0346015,-6.31345e-07,-2.51693e-08,1.16693e-11,70114.7,22.8863], Tmin=(100,'K'), Tmax=(1046.06,'K')), NASAPolynomial(coeffs=[13.7985,0.0140255,-6.73462e-06,1.41395e-09,-1.07535e-13,66271.8,-41.4319], Tmin=(1046.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]=C[C][C]=O(17992)',
    structure = SMILES('[CH2][C]=C[C]=[C][O]'),
    E0 = (785.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.25849,'amu*angstrom^2'), symmetry=1, barrier=(51.9272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25741,'amu*angstrom^2'), symmetry=1, barrier=(51.9023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81157,0.0449805,-4.04526e-05,1.36439e-08,6.06583e-13,94519.9,23.6856], Tmin=(100,'K'), Tmax=(867.791,'K')), NASAPolynomial(coeffs=[11.5058,0.0118814,-3.26548e-06,4.59873e-10,-2.73014e-14,92401.2,-24.2221], Tmin=(867.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH][C]=C[CH][C]=O(17993)',
    structure = SMILES('[CH][C]=C[CH][C]=O'),
    E0 = (801.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,408.972,408.974,408.977,408.98],'cm^-1')),
        HinderedRotor(inertia=(0.411014,'amu*angstrom^2'), symmetry=1, barrier=(48.7818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411002,'amu*angstrom^2'), symmetry=1, barrier=(48.7818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410993,'amu*angstrom^2'), symmetry=1, barrier=(48.7818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91322,0.0370785,-8.66247e-06,-1.3098e-08,6.47325e-12,96474.8,23.7561], Tmin=(100,'K'), Tmax=(1116.78,'K')), NASAPolynomial(coeffs=[11.9087,0.0194327,-9.34659e-06,1.86722e-09,-1.35504e-13,93110.1,-30.6391], Tmin=(1116.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(801.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Cds_S) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC1[C]=O(17994)',
    structure = SMILES('C=C1[CH]C1[C]=O'),
    E0 = (432.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25513,0.0299459,1.5203e-06,-2.40392e-08,1.14425e-11,52035.1,20.8698], Tmin=(100,'K'), Tmax=(985.538,'K')), NASAPolynomial(coeffs=[10.8668,0.0147549,-5.43568e-06,1.01165e-09,-7.31651e-14,49378,-25.4201], Tmin=(985.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C=C=C[C]=O(17995)',
    structure = SMILES('[CH2]C=[C]C=C=O'),
    E0 = (309.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12658,0.043753,-4.93798e-05,3.76275e-08,-1.22826e-11,37282,19.9963], Tmin=(100,'K'), Tmax=(833.058,'K')), NASAPolynomial(coeffs=[5.18955,0.0256982,-1.08424e-05,1.96348e-09,-1.32237e-13,36887.8,6.47597], Tmin=(833.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=C=C[C]=O(17996)',
    structure = SMILES('[CH2]C#C[CH][C]=O'),
    E0 = (503.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2100,2250,500,550,1855,455,950,593.68],'cm^-1')),
        HinderedRotor(inertia=(0.000480806,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142501,'amu*angstrom^2'), symmetry=1, barrier=(35.7339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143459,'amu*angstrom^2'), symmetry=1, barrier=(35.7388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37654,0.0232414,2.2253e-05,-4.84189e-08,2.04699e-11,60573.1,22.0427], Tmin=(100,'K'), Tmax=(982.814,'K')), NASAPolynomial(coeffs=[13.1096,0.00997406,-3.91919e-06,8.2293e-10,-6.56418e-14,56994.5,-37.0253], Tmin=(982.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(503.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CtHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCJC=O) + radical(Propargyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=[C][CH][C]=O(17997)',
    structure = SMILES('[CH2][C]=[C]C=[C][O]'),
    E0 = (785.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.25849,'amu*angstrom^2'), symmetry=1, barrier=(51.9272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.25741,'amu*angstrom^2'), symmetry=1, barrier=(51.9023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81157,0.0449805,-4.04526e-05,1.36439e-08,6.06583e-13,94519.9,23.6856], Tmin=(100,'K'), Tmax=(867.791,'K')), NASAPolynomial(coeffs=[11.5058,0.0118814,-3.26548e-06,4.59873e-10,-2.73014e-14,92401.2,-24.2221], Tmin=(867.791,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C]C[C]=O(17165)',
    structure = SMILES('[CH2][C]=[C]C[C]=O'),
    E0 = (657.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1670,1700,300,440,1855,455,950,371.987],'cm^-1')),
        HinderedRotor(inertia=(0.00121861,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190537,'amu*angstrom^2'), symmetry=1, barrier=(18.701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190421,'amu*angstrom^2'), symmetry=1, barrier=(18.7013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29977,0.0405298,-3.06304e-05,1.10523e-08,-1.64354e-12,79194.6,21.6119], Tmin=(100,'K'), Tmax=(1497.69,'K')), NASAPolynomial(coeffs=[9.81525,0.0204579,-1.05278e-05,2.10412e-09,-1.49905e-13,76943.4,-17.6817], Tmin=(1497.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=[C][CH][C]=O(17998)',
    structure = SMILES('[CH2]C=[C]C=[C][O]'),
    E0 = (547.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81951,0.0406988,-1.68885e-05,-1.38622e-08,1.07318e-11,65917.6,23.1818], Tmin=(100,'K'), Tmax=(902.379,'K')), NASAPolynomial(coeffs=[12.7371,0.0123069,-2.94349e-06,4.0015e-10,-2.50693e-14,63132.9,-32.8802], Tmin=(902.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C][CH]C=O(17999)',
    structure = SMILES('[CH2][C]=[C]C=C[O]'),
    E0 = (545.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5791,0.0435018,-1.46016e-05,-2.37305e-08,1.61967e-11,65699.9,21.3905], Tmin=(100,'K'), Tmax=(893.788,'K')), NASAPolynomial(coeffs=[15.7012,0.00765852,-3.613e-07,-1.05605e-10,9.57888e-15,62082.7,-51.2681], Tmin=(893.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=[C][CH][C]=O(18000)',
    structure = SMILES('C[C]=[C][CH][C]=O'),
    E0 = (668.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,869.371],'cm^-1')),
        HinderedRotor(inertia=(0.0337131,'amu*angstrom^2'), symmetry=1, barrier=(17.9292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781275,'amu*angstrom^2'), symmetry=1, barrier=(17.9631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283648,'amu*angstrom^2'), symmetry=1, barrier=(17.9422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96918,0.0405467,-2.88747e-05,9.29416e-09,-1.16895e-12,80491.4,22.9309], Tmin=(100,'K'), Tmax=(1854.06,'K')), NASAPolynomial(coeffs=[14.7366,0.013002,-6.59013e-06,1.28129e-09,-8.85016e-14,75757,-46.5469], Tmin=(1854.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Cds_S) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C=C[C]=O(16163)',
    structure = SMILES('C=[C]C=C[C]=O'),
    E0 = (330.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.23543,'amu*angstrom^2'), symmetry=1, barrier=(28.4051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23355,'amu*angstrom^2'), symmetry=1, barrier=(28.3617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43864,0.0597968,-7.86998e-05,5.54009e-08,-1.56855e-11,39785,19.3972], Tmin=(100,'K'), Tmax=(860.663,'K')), NASAPolynomial(coeffs=[9.97258,0.0201332,-9.56998e-06,1.85128e-09,-1.30193e-13,38316.1,-20.4933], Tmin=(860.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C1=CC=C1[O](18001)',
    structure = SMILES('C=C1[CH]C=C1[O]'),
    E0 = (245.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26723,0.0160275,7.28256e-05,-1.20897e-07,5.22032e-11,29593.2,15.148], Tmin=(100,'K'), Tmax=(909.654,'K')), NASAPolynomial(coeffs=[19.0986,-0.00023417,4.41134e-06,-9.66423e-10,6.17575e-14,24141.7,-77.5931], Tmin=(909.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2][C]=CC#C[O](18002)',
    structure = SMILES('C=[C][CH][C]=C=O'),
    E0 = (543.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,2120,512.5,787.5,194.187,194.192],'cm^-1')),
        HinderedRotor(inertia=(1.57828,'amu*angstrom^2'), symmetry=1, barrier=(42.235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5783,'amu*angstrom^2'), symmetry=1, barrier=(42.235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87342,0.0453999,-4.56675e-05,2.39369e-08,-5.02902e-12,65448,18.4218], Tmin=(100,'K'), Tmax=(1149.41,'K')), NASAPolynomial(coeffs=[10.5067,0.0153555,-6.45891e-06,1.19558e-09,-8.27037e-14,63463.4,-24.4311], Tmin=(1149.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CCJC=C=O) + radical(Cds_S) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH2][C]=C[C]=C[O](18003)',
    structure = SMILES('[CH2][C]=C[C]=C[O]'),
    E0 = (545.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5791,0.0435018,-1.46016e-05,-2.37305e-08,1.61967e-11,65699.9,21.3905], Tmin=(100,'K'), Tmax=(893.788,'K')), NASAPolynomial(coeffs=[15.7012,0.00765852,-3.613e-07,-1.05605e-10,9.57888e-15,62082.7,-51.2681], Tmin=(893.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[C]=[C][O](18004)',
    structure = SMILES('[CH2]C=C[C]=[C][O]'),
    E0 = (547.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81951,0.0406988,-1.68885e-05,-1.38622e-08,1.07318e-11,65917.6,23.1818], Tmin=(100,'K'), Tmax=(902.379,'K')), NASAPolynomial(coeffs=[12.7371,0.0123069,-2.94349e-06,4.0015e-10,-2.50693e-14,63132.9,-32.8802], Tmin=(902.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=C[C]=[C]O(18005)',
    structure = SMILES('[CH2][C]=C[C]=[C]O'),
    E0 = (643.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180],'cm^-1')),
        HinderedRotor(inertia=(1.59725,'amu*angstrom^2'), symmetry=1, barrier=(36.7239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59261,'amu*angstrom^2'), symmetry=1, barrier=(36.6172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59395,'amu*angstrom^2'), symmetry=1, barrier=(36.6481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919946,0.0577652,-6.49632e-05,3.60301e-08,-7.44669e-12,77543.6,25.4652], Tmin=(100,'K'), Tmax=(1376.47,'K')), NASAPolynomial(coeffs=[14.9524,0.00770797,-3.01529e-07,-1.85096e-10,2.03608e-14,74559.6,-43.5236], Tmin=(1376.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = 'C[C]=C[C]=[C][O](18006)',
    structure = SMILES('C[C]=C[C]=[C][O]'),
    E0 = (667.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03884,'amu*angstrom^2'), symmetry=1, barrier=(23.8849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02856,'amu*angstrom^2'), symmetry=1, barrier=(23.6485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51448,0.054131,-6.48465e-05,4.15729e-08,-1.05298e-11,80329.2,23.1913], Tmin=(100,'K'), Tmax=(969.803,'K')), NASAPolynomial(coeffs=[10.7746,0.0159382,-5.77489e-06,9.66503e-10,-6.23465e-14,78533.1,-21.1995], Tmin=(969.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C]C=[C]O(18007)',
    structure = SMILES('[CH2][C]=[C]C=[C]O'),
    E0 = (643.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180],'cm^-1')),
        HinderedRotor(inertia=(1.59725,'amu*angstrom^2'), symmetry=1, barrier=(36.7239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59261,'amu*angstrom^2'), symmetry=1, barrier=(36.6172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59395,'amu*angstrom^2'), symmetry=1, barrier=(36.6481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919946,0.0577652,-6.49632e-05,3.60301e-08,-7.44669e-12,77543.6,25.4652], Tmin=(100,'K'), Tmax=(1376.47,'K')), NASAPolynomial(coeffs=[14.9524,0.00770797,-3.01529e-07,-1.85096e-10,2.03608e-14,74559.6,-43.5236], Tmin=(1376.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C1C=C1[O](18008)',
    structure = SMILES('C=[C]C1[CH]C1=O'),
    E0 = (475.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0054,0.0347915,-3.51346e-06,-2.33142e-08,1.22518e-11,57220.7,18.7134], Tmin=(100,'K'), Tmax=(964.557,'K')), NASAPolynomial(coeffs=[12.2978,0.0138293,-4.69191e-06,8.45846e-10,-6.10242e-14,54224.8,-35.8074], Tmin=(964.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(CCJC=O) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C=C[O](18009)',
    structure = SMILES('C=[C]C=C=C[O]'),
    E0 = (366.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44293,0.0581897,-5.88991e-05,2.83074e-08,-4.99587e-12,44184.9,22.5048], Tmin=(100,'K'), Tmax=(1656.62,'K')), NASAPolynomial(coeffs=[16.9923,0.00540494,5.08161e-07,-2.72965e-10,2.24497e-14,40461.6,-60.3782], Tmin=(1656.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[CH]C=[C][O](18010)',
    structure = SMILES('[CH]=C=C[CH][C]=O'),
    E0 = (524.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,540,610,2055,1855,455,950,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.4158,'amu*angstrom^2'), symmetry=1, barrier=(32.5521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41458,'amu*angstrom^2'), symmetry=1, barrier=(32.5239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00033,0.0332974,-2.21907e-06,-2.52658e-08,1.27188e-11,63104.8,21.5334], Tmin=(100,'K'), Tmax=(997.245,'K')), NASAPolynomial(coeffs=[14.044,0.00974557,-4.02982e-06,8.37378e-10,-6.53537e-14,59471.8,-42.7088], Tmin=(997.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJC=O) + radical(CCCJ=O) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C[C][C]=O(16161)',
    structure = SMILES('C=[C]C[C][C]=O'),
    E0 = (722.916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,299.46,299.461,299.462],'cm^-1')),
        HinderedRotor(inertia=(0.108831,'amu*angstrom^2'), symmetry=1, barrier=(6.92589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108833,'amu*angstrom^2'), symmetry=1, barrier=(6.92593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2753,'amu*angstrom^2'), symmetry=1, barrier=(17.5193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62474,0.056845,-8.49176e-05,7.1764e-08,-2.41607e-11,87027.8,24.0783], Tmin=(100,'K'), Tmax=(832.116,'K')), NASAPolynomial(coeffs=[7.27816,0.02262,-1.05157e-05,1.97524e-09,-1.34907e-13,86331,-0.690763], Tmin=(832.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C[CH][CH][C]=O(14871)',
    structure = SMILES('[CH]C=C[CH][C]=O'),
    E0 = (563.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,385.348,394.161,394.408,395.185],'cm^-1')),
        HinderedRotor(inertia=(0.430355,'amu*angstrom^2'), symmetry=1, barrier=(48.5288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438793,'amu*angstrom^2'), symmetry=1, barrier=(48.4166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.442046,'amu*angstrom^2'), symmetry=1, barrier=(48.5268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95762,0.0324086,1.60477e-05,-4.17637e-08,1.69483e-11,67870.9,23.1189], Tmin=(100,'K'), Tmax=(1037.9,'K')), NASAPolynomial(coeffs=[12.6177,0.0207039,-9.49536e-06,1.91583e-09,-1.42082e-13,64075.7,-36.3294], Tmin=(1037.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(AllylJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C]C[CH][C]=O(16166)',
    structure = SMILES('[CH]=[C]C[CH][C]=O'),
    E0 = (689.327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1685,370,1855,455,950,3120,650,792.5,1650,309.189],'cm^-1')),
        HinderedRotor(inertia=(0.187291,'amu*angstrom^2'), symmetry=1, barrier=(12.7217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187692,'amu*angstrom^2'), symmetry=1, barrier=(12.7303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0113573,'amu*angstrom^2'), symmetry=1, barrier=(74.5251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79282,0.0506005,-6.3105e-05,4.45979e-08,-1.27993e-11,82984.4,23.7372], Tmin=(100,'K'), Tmax=(848.507,'K')), NASAPolynomial(coeffs=[8.37481,0.0195721,-8.25322e-06,1.50146e-09,-1.01672e-13,81867.4,-6.93589], Tmin=(848.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C=C[O](18011)',
    structure = SMILES('[CH][C]=CC=C[O]'),
    E0 = (599.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0564,'amu*angstrom^2'), symmetry=1, barrier=(47.2808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0451,'amu*angstrom^2'), symmetry=1, barrier=(47.0208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37422,0.0458517,-8.79841e-06,-3.03402e-08,1.76326e-11,72159.6,22.0779], Tmin=(100,'K'), Tmax=(930.652,'K')), NASAPolynomial(coeffs=[16.2986,0.0119822,-3.00724e-06,4.68084e-10,-3.37592e-14,68070.6,-55.8954], Tmin=(930.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C=[C]O(18012)',
    structure = SMILES('[CH][C]=CC=[C]O'),
    E0 = (697.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05807,'amu*angstrom^2'), symmetry=1, barrier=(47.3191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04982,'amu*angstrom^2'), symmetry=1, barrier=(47.1294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05005,'amu*angstrom^2'), symmetry=1, barrier=(47.1347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635118,0.0608593,-6.08081e-05,3.01296e-08,-5.65854e-12,84007,26.4526], Tmin=(100,'K'), Tmax=(1475.52,'K')), NASAPolynomial(coeffs=[16.3942,0.0106694,-2.19314e-06,2.15882e-10,-8.9969e-15,80169.5,-52.9513], Tmin=(1475.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(697.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CC[C]=O(17162)',
    structure = SMILES('[CH][C]=CC[C]=O'),
    E0 = (639.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,396.374,396.374,396.374,396.375],'cm^-1')),
        HinderedRotor(inertia=(0.447074,'amu*angstrom^2'), symmetry=1, barrier=(49.8442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447072,'amu*angstrom^2'), symmetry=1, barrier=(49.8442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.447072,'amu*angstrom^2'), symmetry=1, barrier=(49.8442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96663,0.0421579,-2.52349e-05,6.61558e-09,-6.68807e-13,76968.1,23.1997], Tmin=(100,'K'), Tmax=(2287.04,'K')), NASAPolynomial(coeffs=[18.1026,0.0139368,-6.726e-06,1.22039e-09,-7.90624e-14,69587.2,-67.9962], Tmin=(2287.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=C[C][CH2](17838)',
    structure = SMILES('[CH]C=[C][CH2]'),
    E0 = (731.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,296.375,297.489,297.589],'cm^-1')),
        HinderedRotor(inertia=(0.828676,'amu*angstrom^2'), symmetry=1, barrier=(51.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815826,'amu*angstrom^2'), symmetry=1, barrier=(51.0347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,88067.9,15.929], Tmin=(100,'K'), Tmax=(2065.82,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,84627.8,-30.2813], Tmin=(2065.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(731.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]C1[CH]C1=O(18013)',
    structure = SMILES('[CH2][C]C1[CH]C1=O'),
    E0 = (764.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77827,0.0380689,-4.00048e-06,-2.86265e-08,1.55433e-11,91985.9,20.9409], Tmin=(100,'K'), Tmax=(946.365,'K')), NASAPolynomial(coeffs=[14.5677,0.0107333,-3.02735e-06,5.2434e-10,-3.92864e-14,88368.6,-46.3779], Tmin=(946.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(764.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CCJC=O) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'O=[C]C1[CH][C]C1(18014)',
    structure = SMILES('O=[C]C1[CH][C]C1'),
    E0 = (708.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49483,0.0194653,3.86395e-05,-6.55623e-08,2.66692e-11,85260,22.2713], Tmin=(100,'K'), Tmax=(963.688,'K')), NASAPolynomial(coeffs=[12.2321,0.0130801,-4.39249e-06,8.51087e-10,-6.60479e-14,81803,-32.5442], Tmin=(963.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCC=O) + radical(CCJ2_triplet) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C][CH]C1[C]C1(18015)',
    structure = SMILES('O=[C][CH]C1[C]C1'),
    E0 = (678.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11194,0.0304427,1.05263e-05,-4.09773e-08,1.95883e-11,81668,20.4157], Tmin=(100,'K'), Tmax=(935.836,'K')), NASAPolynomial(coeffs=[13.3253,0.0105767,-2.61196e-06,4.24989e-10,-3.19114e-14,78340.4,-39.5044], Tmin=(935.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(678.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = '[CH][C][CH2](17602)',
    structure = SMILES('[CH][C][CH2]'),
    E0 = (981.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,917.882,918.218,918.236],'cm^-1')),
        HinderedRotor(inertia=(0.00702614,'amu*angstrom^2'), symmetry=1, barrier=(4.20306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182792,'amu*angstrom^2'), symmetry=1, barrier=(4.20276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96474,0.0202917,-1.68078e-05,7.2211e-09,-1.21436e-12,118107,14.5039], Tmin=(100,'K'), Tmax=(1441.81,'K')), NASAPolynomial(coeffs=[8.09928,0.00604662,-1.9874e-06,3.68258e-10,-2.60955e-14,116627,-12.1458], Tmin=(1441.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(981.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH][C]=C[O](18016)',
    structure = SMILES('[CH2][C][CH][C]=C[O]'),
    E0 = (924.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,421.343,421.401,421.415,421.468],'cm^-1')),
        HinderedRotor(inertia=(0.291513,'amu*angstrom^2'), symmetry=1, barrier=(36.7435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000949007,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291606,'amu*angstrom^2'), symmetry=1, barrier=(36.7402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46024,0.0478458,-3.29672e-05,1.31098e-09,4.87864e-12,111324,24.4596], Tmin=(100,'K'), Tmax=(965.036,'K')), NASAPolynomial(coeffs=[14.8992,0.0104924,-3.42931e-06,6.09373e-10,-4.40196e-14,107876,-44.3244], Tmin=(965.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(924.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJ) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]C1C=[C]O1(18017)',
    structure = SMILES('[CH2][C]C1C=[C]O1'),
    E0 = (831.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80868,0.0298333,3.20445e-05,-7.58643e-08,3.46283e-11,100154,22.7075], Tmin=(100,'K'), Tmax=(933.977,'K')), NASAPolynomial(coeffs=[19.2483,0.00201963,1.42974e-06,-2.74023e-10,1.07316e-14,94851.5,-71.1827], Tmin=(933.977,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(RCCJ) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = 'O=C1[CH][CH][C]C1(18018)',
    structure = SMILES('[O]C1=C[CH][C]C1'),
    E0 = (520.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4698,0.0157898,5.93278e-05,-9.4301e-08,3.88427e-11,62704.9,16.5144], Tmin=(100,'K'), Tmax=(939.409,'K')), NASAPolynomial(coeffs=[14.9072,0.00829204,-1.28944e-06,2.31202e-10,-2.39124e-14,58362.2,-53.3881], Tmin=(939.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=C(C)OJ) + radical(cyclopentene-allyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]1[CH]C=[C]OC1(18019)',
    structure = SMILES('[C]1[CH]C=[C]OC1'),
    E0 = (674.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6564,0.0128097,6.4746e-05,-9.93316e-08,4.12643e-11,81209.9,14.9439], Tmin=(100,'K'), Tmax=(915.111,'K')), NASAPolynomial(coeffs=[13.3786,0.0092819,-5.10882e-07,-3.85165e-11,-8.78075e-17,77432.8,-45.7482], Tmin=(915.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(C=CJO)"""),
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
    label = '[CH2]C=[C][CH][C][O](18020)',
    structure = SMILES('[CH2][CH][C][CH][C]=O'),
    E0 = (924.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,405.974,406.468,2666.75],'cm^-1')),
        HinderedRotor(inertia=(0.0662472,'amu*angstrom^2'), symmetry=1, barrier=(7.75709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472488,'amu*angstrom^2'), symmetry=1, barrier=(55.2109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109486,'amu*angstrom^2'), symmetry=1, barrier=(55.2264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109476,'amu*angstrom^2'), symmetry=1, barrier=(55.2208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74778,0.0546482,-8.23344e-05,7.07153e-08,-2.3697e-11,111270,27.6646], Tmin=(100,'K'), Tmax=(882.94,'K')), NASAPolynomial(coeffs=[6.18541,0.0233417,-1.0117e-05,1.81739e-09,-1.20105e-13,110923,9.28087], Tmin=(882.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(924.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C[CH][CH][C][O](18021)',
    structure = SMILES('[CH][CH][CH]C=[C][O]'),
    E0 = (915.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,868.487,868.491,868.491,868.498,868.498,868.562],'cm^-1')),
        HinderedRotor(inertia=(0.0057093,'amu*angstrom^2'), symmetry=1, barrier=(3.05616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132966,'amu*angstrom^2'), symmetry=1, barrier=(3.05715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132894,'amu*angstrom^2'), symmetry=1, barrier=(3.05551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03737,0.0424217,-3.67534e-05,1.74581e-08,-3.44635e-12,110227,26.3987], Tmin=(100,'K'), Tmax=(1193.01,'K')), NASAPolynomial(coeffs=[8.80961,0.0197153,-8.20405e-06,1.50436e-09,-1.03184e-13,108611,-7.46861], Tmin=(1193.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(915.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJC) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = 'C=C1[CH]C1[C][O](18022)',
    structure = SMILES('C=C1[CH]C1[C][O]'),
    E0 = (823.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94774,0.0381306,-1.8867e-05,-2.4354e-09,3.33473e-12,99084.4,18.0417], Tmin=(100,'K'), Tmax=(1082.6,'K')), NASAPolynomial(coeffs=[11.2972,0.0165743,-6.99534e-06,1.33591e-09,-9.52581e-14,96299,-31.3211], Tmin=(1082.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(823.171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CCOJ) + radical(Allyl_S) + radical(CH2_triplet)"""),
)

species(
    label = 'C=C1[CH][CH][C]O1(18023)',
    structure = SMILES('[CH2]C1=C[CH][C]O1'),
    E0 = (566.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7716,0.0309622,3.31831e-05,-8.27604e-08,3.9694e-11,68232.2,12.9654], Tmin=(100,'K'), Tmax=(896.969,'K')), NASAPolynomial(coeffs=[19.7505,-0.000576502,4.58853e-06,-1.05471e-09,7.20169e-14,63050.4,-82.7243], Tmin=(896.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(2,3-Dihydrofuran) + radical(C=CCJCO) + radical(C=C(O)CJ) + radical(CH2_triplet)"""),
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
    label = 'C=[C][CH]C1[C]O1(18024)',
    structure = SMILES('C=[C][CH]C1[C]O1'),
    E0 = (792.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.671372,0.0596522,-6.48209e-05,3.49243e-08,-6.94644e-12,95451.8,20.066], Tmin=(100,'K'), Tmax=(1468.81,'K')), NASAPolynomial(coeffs=[14.9828,0.00830552,2.51132e-07,-3.45545e-10,3.28012e-14,92582.3,-49.9375], Tmin=(1468.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(792.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CCJCO) + radical(CH2_triplet) + radical(Cds_S)"""),
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
    E0 = (582.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1094.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (996.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1013.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (585.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (605.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (722.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (780.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (996.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (840.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (763.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (735.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (822.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (582.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (590.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (770.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (802.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (731.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (775.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (836.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (836.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (792.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (585.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (605.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (750.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (905.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (777.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (834.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (632.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (848.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (848.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1225.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (764.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (726.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (678.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (632.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (1202.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (947.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (831.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (698.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (674.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (1279.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (947.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (938.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (823.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (656.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (1227.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (795.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['HCCO(2227)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2][C]=C[C][C]=O(17992)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH][C]=C[CH][C]=O(17993)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2]C1=CC1[C]=O(17994)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2]C=C=C[C]=O(17995)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C]=C=C[C]=O(17996)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C][O](6861)', 'C3H3(5450)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(72.9547,'m^3/(mol*s)'), n=1.66457, Ea=(16.7701,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;YJ] for rate rule [Ct_Ct;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]=[C][CH][C]=O(17997)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=[C]C[C]=O(17165)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.28974e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2]C=[C][CH][C]=O(17998)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.97782e+08,'s^-1'), n=1.48417, Ea=(181.505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2][C]=[C][CH]C=O(17999)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[C]=[C][CH][C]=O(18000)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.09178e+08,'s^-1'), n=1.56859, Ea=(153.918,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['C=[C]C=C[C]=O(16163)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2]C1=CC=C1[O](18001)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2][C]=CC#C[O](18002)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HCCO(2227)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2][C]=C[C]=C[O](18003)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.44274e+08,'s^-1'), n=1.26608, Ea=(149.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2]C=C[C]=[C][O](18004)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C]=C[C]=[C]O(18005)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[C]=C[C]=[C][O](18006)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(115.297,'s^-1'), n=2.99825, Ea=(169.04,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]=[C]C=[C]O(18007)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1194.19,'s^-1'), n=2.56519, Ea=(148.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;O_H_out] + [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['C=[C]C1C=C1[O](18008)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['C=[C]C=C=C[O](18009)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', 'C#C[CH]C=[C][O](18010)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C[C][C]=O(16161)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.28974e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH]=C[CH][CH][C]=O(14871)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]C[CH][C]=O(16166)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C][CH]C=C[O](18011)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C][CH]C=[C]O(18012)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH][C]=CC[C]=O(17162)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]=O(2355)', '[CH]=C[C][CH2](17838)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2][C]C1[CH]C1=O(18013)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(182.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['O=[C]C1[CH][C]C1(18014)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.14066e+09,'s^-1'), n=0.807397, Ea=(144.515,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['O=[C][CH]C1[C]C1(18015)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(96.1104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 95.6 to 96.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CO(2039)', '[CH]=C[C][CH2](17838)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction37',
    reactants = ['HCCO(2227)', '[CH][C][CH2](17602)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][C][CH][C]=C[O](18016)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[CH2][C]C1C=[C]O1(18017)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.73504e+09,'s^-1'), n=0.685238, Ea=(249.663,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 249.4 to 249.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['O=C1[CH][CH][C]C1(18018)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.49459e+10,'s^-1'), n=0.314867, Ea=(116.337,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['[C]1[CH]C=[C]OC1(18019)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(92.4196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;multiplebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 89.7 to 92.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH][C][O](10218)', 'C3H3(5450)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C=[C][CH][C][O](18020)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C[CH][CH][C][O](18021)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['C=C1[CH]C1[C][O](18022)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(241.233,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['C=C1[CH][CH][C]O1(18023)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['H2CC(T)(1341)', '[CH]=C[C][O](10430)'],
    products = ['C=[C][CH][CH][C]=O(16165)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C][CH][CH][C]=O(16165)'],
    products = ['C=[C][CH]C1[C]O1(18024)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.30963e+11,'s^-1'), n=0.419784, Ea=(213.014,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R4;multiplebond_intra;radadd_intra_O] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = 'PDepNetwork #4052',
    isomers = [
        'C=[C][CH][CH][C]=O(16165)',
    ],
    reactants = [
        ('HCCO(2227)', 'C3H3(5450)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4052',
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


import math
import numpy as np

def calcola_MPP_Ac_Power(dayNumber, localStandardTimeInHoursFromMidnight, DNI, gDiffH, windSpeed, Ta, SITE, PV, SUBARRAY, INV):
    # Calcola le grandezze astronomiche
    delta, omega, I0, thetaZeta, m = calcola_grandezze_astronomiche(dayNumber, localStandardTimeInHoursFromMidnight, SITE['GMT'], SITE['latitude'], SITE['longitude'], SITE['z'])
    
    # Calcola l'angolo di incidenza su piano inclinato
    theta = calcola_angolo_incidenza_su_piano_inclinato(delta, omega, SITE['latitude'], SUBARRAY['betaTilt'], SUBARRAY['gammaTilt'])
    
    # Calcola l'irraggiamento diffuso inclinato con il modello di Perez
    Ediff_tilted = calcola_Irr_Tilted_Con_Perez(DNI, gDiffH, I0, thetaZeta, theta, SUBARRAY['betaTilt'], m)
    
    # Calcola la tensione e la corrente al punto di massima potenza del subarray
    Vdc, Idc = calcola_MPP_subArray(DNI, Ediff_tilted, theta, m, windSpeed, Ta, SUBARRAY['Ns'], SUBARRAY['Np'], PV)
    
    # Calcola la potenza continua
    Pdc = Vdc * Idc
    
    # Calcola la potenza AC
    AcPow = SUBARRAY['TotDerate'] * calcola_Pac_inverter(Pdc, Vdc, INV)
    
    return AcPow


def calcola_Pac_inverter(Pdc, Vdc, INV):
    """
    Calcola la potenza AC in uscita dall'inverter.

    Parametri:
    Pdc (float): Potenza continua in ingresso all'inverter [W]
    Vdc (float): Tensione continua in ingresso all'inverter [V]
    INV (dict): Dizionario con i parametri dell'inverter.

    Returns:
    Pac (float): Potenza AC in uscita dall'inverter [W]
    """
    Pac = 0
    if Pdc > 0:
        A = INV['Pdc0'] * (1 + INV['C1'] * (Vdc - INV['Vdc0']))
        B = INV['Ps0'] * (1 + INV['C2'] * (Vdc - INV['Vdc0']))
        C = INV['C0'] * (1 + INV['C3'] * (Vdc - INV['Vdc0']))
        Pac = max(0, (INV['Pac0'] / (A - B) - C * (A - B)) * (Pdc - B) + C * (Pdc - B) ** 2)
    return Pac

def calcola_MPP_subArray(DNI, Ediff_tilted, AOI, AMa, WS, Ta, Ns, Np, PV):
    """
    Calcola la tensione e la corrente al punto di massima potenza del subarray.

    Parametri:
    DNI (float): Irraggiamento diretto normale [W/m^2]
    Ediff_tilted (float): Irraggiamento diffuso inclinato [W/m^2]
    AOI (float): Angolo di incidenza su piano inclinato [gradi]
    AMa (float): Air Mass (AM) airmass [unitless]
    WS (float): Velocità del vento [m/s]
    Ta (float): Temperatura ambiente [°C]
    Ns (int): Numero di moduli in serie in una singola stringa [-]
    Np (int): Numero di stringhe in parallelo in un subarray [-]
    PV (dict): Dizionario con i parametri del modulo fotovoltaico.

    Returns:
    Vmp (float): Tensione al punto di massima potenza del subarray [V]
    Imp (float): Corrente al punto di massima potenza del subarray [A]
    """
    VmpMod, ImpMod = calcola_MPP_modulo(PV, DNI, Ediff_tilted, AOI, AMa, WS, Ta)
    Vmp = Ns * VmpMod
    Imp = Np * ImpMod
    return Vmp, Imp

def calcola_MPP_modulo(PV, DNI, Ediff_tilted, AOI, AMa, WS, Ta):
    """
    Calcola la tensione e la corrente al punto di massima potenza del modulo fotovoltaico.

    Parametri:
    PV (dict): Dizionario con i parametri del modulo fotovoltaico.
    DNI (float): Irraggiamento diretto normale [W/m^2]
    Ediff_tilted (float): Irraggiamento diffuso inclinato [W/m^2]
    AOI (float): Angolo di incidenza su piano inclinato [gradi]
    AMa (float): Air Mass (AM) airmass [-]
    WS (float): Velocità del vento [m/s]
    Ta (float): Temperatura ambiente [°C]

    Returns:
    Vmp (float): Tensione al punto di massima potenza del modulo [V]
    Imp (float): Corrente al punto di massima potenza del modulo [A]
    """
    k = 1.38066e-23  # Boltzmann [J/K]
    q = 1.60218e-19  # Electron charge [coulomb]

    Vmp = 0
    Imp = 0
    Ebeam_tilted = DNI * math.cos(AOI * math.pi / 180)
    f1 = PV['a0'] + PV['a1'] * AMa + PV['a2'] * AMa ** 2 + PV['a3'] * AMa ** 3 + PV['a4'] * AMa ** 4  # spectral correction [-]
    f2 = PV['b0'] + PV['b1'] * AOI + PV['b2'] * AOI ** 2 + PV['b3'] * AOI ** 3 + PV['b4'] * AOI ** 4 + PV['b5'] * AOI ** 5  # Incidence angle modifier [-]
    Ee = f1 * (f2 * Ebeam_tilted + PV['fd'] * Ediff_tilted) / PV['E0']  # Effective solar irradiance [-]
    
    if Ee > 0:
        E = Ebeam_tilted + Ediff_tilted  # Solar irradiance on module plane [W/m^2]
        Tm = E * math.exp(PV['a'] + PV['b'] * WS) + Ta  # Module back surface temperature [°C]
        Tc = Tm + E / PV['E0'] * PV['DTC']  # Cells temperature inside the module [°C]
        deltaT = PV['n'] * k * (Tc + 273.15) / q  # Cell's thermal voltage [V]
        betaVmp = PV['betaVmp0'] + PV['m_betaVmp'] * (1 - Ee)  # Temp. coefficient for MPP Voltage [V/°C]
        Vmp = PV['Vmp0'] + PV['C2'] * PV['Ns'] * deltaT * math.log(Ee) + PV['C3'] * PV['Ns'] * (deltaT * math.log(Ee)) ** 2 + betaVmp * (Tc - PV['T0'])  # MPP Voltage [V]
        Imp = PV['Imp0'] * (PV['C0'] * Ee + PV['C1'] * Ee ** 2) * (1 + PV['alphaImp'] * (Tc - PV['T0']))
    
    return Vmp, Imp

def calcola_Irr_Tilted_Con_Perez(DNI, gDiffH, I0, thetaZeta, theta, betaTilt, m):
    """
    Calcola l'irraggiamento diffuso inclinato usando il modello di Perez.

    Parametri:
    DNI (float): Irraggiamento diretto normale [W/m^2]
    gDiffH (float): Irraggiamento diffuso orizzontale [W/m^2]
    I0 (float): Irraggiamento extraterrestre [W/m^2]
    thetaZeta (float): Angolo tra il sole e la perpendicolare al piano inclinato [gradi]
    theta (float): Angolo di incidenza del sole sul piano inclinato [gradi]
    betaTilt (float): Angolo di inclinazione del piano inclinato [gradi]
    m (float): Coefficiente m del modello di Perez

    Returns:
    Ediff_tilted (float): Irraggiamento diffuso inclinato [W/m^2]
    """
    Ediff_tilted = 0
    if 0 < thetaZeta < 90 and gDiffH > 0:
        DELTA = m * gDiffH / I0  # Sky's brightness, eq. 2
        thetaZeta = math.radians(thetaZeta)
        skyClearness = ((gDiffH + DNI) / gDiffH + 1.041 * thetaZeta**3) / (1 + 1.041 * thetaZeta**3)  # epsilon (sky clearness), eq 1
        eBin = computeSkyClearnessLevel(skyClearness)  # Discrete sky clearness, TAB 1
        perezModelCoeff = computePerezModelCoeff(eBin)  # TAB 6
        F1 = perezModelCoeff[0] + perezModelCoeff[1] * DELTA + perezModelCoeff[2] * thetaZeta  # Circumsolar Brightening Coefficient , TAB 6
        F2 = perezModelCoeff[3] + perezModelCoeff[4] * DELTA + perezModelCoeff[5] * thetaZeta  # Horizon Brightening Coefficient , TAB 6
        theta = math.radians(theta)
        a = max(0, math.cos(theta))  # pag 281
        b = max(0.087, math.cos(thetaZeta))  # pag 281
        betaTilt = math.radians(betaTilt)
        f = (1 - F1) * (1 + math.cos(betaTilt)) / 2 + F1 * a / b + F2 * math.sin(betaTilt)  # (eq. 9) , pag 281
        Ediff_tilted = f * gDiffH
    
    return Ediff_tilted

def computePerezModelCoeff(eBin):
    """
    Calcola i coefficienti del modello di Perez in base al valore di eBin.

    Parametri:
    eBin (int): Valore di eBin [1-8]

    Returns:
    array (list): Lista dei coefficienti del modello di Perez [6 coefficienti]
    """
    array = [0] * 6

    if eBin == 1:
        array[0] = -0.008
        array[1] = 0.588
        array[2] = -0.062
        array[3] = -0.060
        array[4] = 0.072
        array[5] = -0.022
    elif eBin == 2:
        array[0] = 0.130
        array[1] = 0.683
        array[2] = -0.151
        array[3] = -0.019
        array[4] = 0.066
        array[5] = -0.029
    elif eBin == 3:
        array[0] = 0.330
        array[1] = 0.487
        array[2] = -0.221
        array[3] = 0.055
        array[4] = -0.064
        array[5] = -0.026
    elif eBin == 4:
        array[0] = 0.568
        array[1] = 0.187
        array[2] = -0.295
        array[3] = 0.109
        array[4] = -0.152
        array[5] = -0.014
    elif eBin == 5:
        array[0] = 0.873
        array[1] = -0.392
        array[2] = -0.362
        array[3] = 0.226
        array[4] = -0.462
        array[5] = 0.001
    elif eBin == 6:
        array[0] = 1.132
        array[1] = -1.237
        array[2] = -0.412
        array[3] = 0.288
        array[4] = -0.823
        array[5] = 0.056
    elif eBin == 7:
        array[0] = 1.060
        array[1] = -1.600
        array[2] = -0.359
        array[3] = 0.264
        array[4] = -1.127
        array[5] = 0.131
    elif eBin == 8:
        array[0] = 0.678
        array[1] = -0.327
        array[2] = -0.250
        array[3] = 0.156
        array[4] = -1.377
        array[5] = 0.251

    return array

def computeSkyClearnessLevel(skyClearness):
    """
    Calcola il livello di chiarezza del cielo (eBin) in base al valore di skyClearness.

    Parametri:
    skyClearness (float): Valore di chiarezza del cielo

    Returns:
    eBin (int): Livello di chiarezza del cielo (eBin)
    """
    if skyClearness < 1.065:
        eBin = 1
    elif 1.065 <= skyClearness < 1.230:
        eBin = 2
    elif 1.230 <= skyClearness < 1.500:
        eBin = 3
    elif 1.500 <= skyClearness < 1.950:
        eBin = 4
    elif 1.950 <= skyClearness < 2.800:
        eBin = 5
    elif 2.800 <= skyClearness < 4.500:
        eBin = 6
    elif 4.500 <= skyClearness < 6.200:
        eBin = 7
    elif skyClearness >= 6.200:
        eBin = 8
    else:
        eBin = 0

    return eBin

def calcola_grandezze_astronomiche(d, localStandardTimeInHoursFromMidnight, GMT, latitude, longitude, z):
    """
    Calcola varie grandezze astronomiche come la declinazione del Sole, l'angolo orario del Sole, la radiazione solare e altro.

    Parametri:
    d (float): Numero del giorno dell'anno
    localStandardTimeInHoursFromMidnight (float): Ore dal mezzanotte del tempo standard locale
    GMT (float): Fuso orario (Meridiano di Greenwich)
    latitude (float): Latitudine del luogo in gradi
    longitude (float): Longitudine del luogo in gradi
    z (float): Altitudine sopra il livello del mare

    Returns:
    delta (float): Declinazione del Sole in gradi
    omega (float): Angolo orario del Sole in gradi
    I0 (float): Radiazione solare extraterrestre orizzontale in W/m^2
    thetaZeta (float): Angolo zenitale del Sole in gradi
    m (float): Massa d'aria ottica del Sole
    """
    gamma = 2 * math.pi * (d - 1.0) / 365.0  # Angolo giornaliero [rad]
    delta = 6.918e-3 - 3.99912e-1 * math.cos(gamma) + 7.0257e-2 * math.sin(gamma) - 6.758e-3 * math.cos(2 * gamma) + 9.07e-4 * math.sin(2 * gamma) - 2.697e-3 * math.cos(3 * gamma) + 1.48e-3 * math.sin(3 * gamma)  # Declinazione del Sole [rad]
    delta = delta * 180 / math.pi

    EQUATION_OF_TIME = 229.18 * (7.5e-5 + 1.868e-3 * math.cos(gamma) - 3.2077e-2 * math.sin(gamma) - 1.4615e-2 * math.cos(2 * gamma) - 4.089e-2 * math.sin(2 * gamma))
    #localStandardTimeInMinsFromMidnight = localStandardTimeInHoursFromMidnight * 60
    localStandardTimeInMinsFromMidnight = localStandardTimeInHoursFromMidnight.hour * 60 + localStandardTimeInHoursFromMidnight.minute
    stdLongitude = GMT * 15 * math.pi / 180
    locLongitude = longitude * math.pi / 180
    tlat = localStandardTimeInMinsFromMidnight + 24 * 60 / (2 * math.pi) * (stdLongitude - locLongitude) + EQUATION_OF_TIME  # Ora solare locale
    omega = (12 * 60 - tlat) * 2 * math.pi / (24 * 60)  # Angolo orario [rad]
    omega = omega * 180 / math.pi

    thetaZeta = calcola_angolo_incidenza_su_piano_inclinato(delta, omega, latitude, 0, 0)  # Angolo zenitale [deg]
    p_p0_ratio = math.exp(-z / 8434.5)

    m = 0
    I0 = 0
    if thetaZeta < 90:
        m = 1 / (math.cos(thetaZeta * math.pi / 180) + 0.50572 / ((96.07995 - thetaZeta) ** 1.6364)) * p_p0_ratio  # Massa d'aria ottica
        ECCENTRICITY_CORRECTION = (1.00011 + 3.4221e-2 * math.cos(gamma) + 1.28e-3 * math.sin(gamma) + 7.19e-4 * math.cos(2 * gamma) + 7.7e-5 * math.sin(2 * gamma))
        I0 = 1367.0 * ECCENTRICITY_CORRECTION * math.cos(thetaZeta * math.pi / 180)  # Radiazione solare extraterrestre

    return delta, omega, I0, thetaZeta, m

def calcola_angolo_incidenza_su_piano_inclinato(delta, omega, latitude, betaTilt, gammaTilt):
    """
    Calcola l'angolo di incidenza del sole su un piano inclinato.

    Parametri:
    delta (float): Declinazione del Sole in gradi
    omega (float): Angolo orario del Sole in gradi
    latitude (float): Latitudine del luogo in gradi
    betaTilt (float): Angolo di inclinazione del piano in gradi
    gammaTilt (float): Angolo di azimuth del piano in gradi

    Returns:
    theta (float): Angolo zenitale rispetto al piano inclinato in gradi
    """
    delta = delta * math.pi / 180
    omega = omega * math.pi / 180
    latitude = latitude * math.pi / 180
    betaTilt = betaTilt * math.pi / 180
    gammaTilt = 180 - gammaTilt
    gammaTilt = gammaTilt * math.pi / 180

    A_COS_THETA = (math.sin(latitude) * math.cos(betaTilt) - math.cos(latitude) * math.sin(betaTilt) * math.cos(gammaTilt)) * math.sin(delta)
    B_COS_THETA = ((math.cos(latitude) * math.cos(betaTilt) + math.sin(latitude) * math.sin(betaTilt) * math.cos(gammaTilt)) * math.cos(delta))
    C_COS_THETA = (math.cos(delta) * math.sin(betaTilt) * math.sin(gammaTilt))
    theta = math.acos(A_COS_THETA + B_COS_THETA * math.cos(omega) + C_COS_THETA * math.sin(omega))  # [rad]
    theta = theta * 180 / math.pi

    return theta

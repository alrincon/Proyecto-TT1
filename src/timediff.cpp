#include "../include/timediff.h"


/**
 * @brief Calculates various time differences between standard time systems.
 *
 * This function computes the differences between several time standards based on
 * the input differences between UT1 and UTC, and between TAI and UTC. The time standards
 * involved include Universal Time 1 (UT1), Coordinated Universal Time (UTC), International
 * Atomic Time (TAI), Global Positioning System time (GPS), and Terrestrial Time (TT).
 *
 * @param UT1_UTC Difference between UT1 and UTC (input).
 * @param TAI_UTC Difference between TAI and UTC (input).
 * @param UT1_TAI Output parameter for the difference between UT1 and TAI.
 * @param UTC_GPS Output parameter for the difference between UTC and GPS.
 * @param UT1_GPS Output parameter for the difference between UT1 and GPS.
 * @param TT_UTC Output parameter for the difference between TT (Terrestrial Time) and UTC.
 * @param GPS_UTC Output parameter for the difference between GPS and UTC.
 */
void  timediff (double UT1_UTC, double TAI_UTC, double& UT1_TAI, double& UTC_GPS,  double& UT1_GPS,  double& TT_UTC,  double& GPS_UTC){
    double TT_TAI  = +32.184;          // TT-TAI time difference [s]

    double GPS_TAI = -19.0;            // GPS-TAI time difference [s]

    double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]

    UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

    UT1_GPS = UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]

    TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

    GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]
}
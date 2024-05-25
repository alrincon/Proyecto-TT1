#ifndef PROYECTO_TIMEDIFF_H
#define PROYECTO_TIMEDIFF_H

#include <iostream>

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
void  timediff (double UT1_UTC, double TAI_UTC, double& UT1_TAI, double& UTC_GPS,  double& UT1_GPS,  double& TT_UTC,  double& GPS_UTC);


#endif //PROYECTO_TIMEDIFF_H

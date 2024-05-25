#ifndef SIGN__H
#define SIGN__H

/**
 * @brief Returns the absolute value of `a` with the sign of `b`.
 *
 * This function determines the absolute value of the first argument `a` and
 * returns it with the sign that corresponds to the sign of the second argument `b`.
 *
 * @param a The value whose absolute value is considered.
 * @param b The value whose sign is used to determine the sign of the result.
 * @return double The absolute value of `a` with the sign of `b`.
 *        If `b` is zero or positive, the result is positive.
 *        If `b` is negative, the result is negative.
 */
double sign(double a, double b);

#endif

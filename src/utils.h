#ifndef ROADEF2020_UTILS_H
#define ROADEF2020_UTILS_H

#include <cstddef>
#include <cmath>
#include <string>
#include <tuple>
#include <optional>
#include <type_traits>


namespace orcs::utils {

    /**
     * Compare two numbers. This function considers a threshold value to compare the numbers. If the difference
     * between they are less then threshold value, then the numbers are considered as equals.
     * @param first The first number.
     * @param second The second number.
     * @param threshold Threshold value used when compared numbers for equality. If not set, then 1E-5 is used.
     * @return  Return -1 if the first number is less then the second number, 0 if both numbers are equal, 1 if
     * the first number is greater than the second number.
     */
    inline int compare(double first, double second, double threshold = 1e-5);

    /**
     * String formatter.
     * @param str A string specifying who to interpret the data.
     * @param args Arguments specifying data to print.
     * @return A formatted string.
     */
    template <class... T>
    std::string format(const std::string& str, T... args);

    /**
     * Check if two intervals, x = [x1, x2] and y = [y1, y2], overlaps.
     * @param x1 Lower bound of interval x.
     * @param x2 Upper bound of interval x.
     * @param y1 Lower bound of interval y.
     * @param y1 Upper bound of interval y.
     * @return True if the intervals overpas, false otherwise.
     */
    template <class T>
    bool overlap(const T& x1, const T& x2, const T& y1, const T& y2);

    /**
     * Get the intersection between two interval x = [x1, x2] and y = [y1, y2].
     * @param x1 Lower bound of interval x.
     * @param x2 Upper bound of interval x.
     * @param y1 Lower bound of interval y.
     * @param y1 Upper bound of interval y.
     * @return An optional return with a tuple containing the range in which the two interval intersect, otherwise, the
     * return is empty.
     */
    template <class T>
    std::optional< std::tuple<T, T> > intersection(const T& x1, const T& x2, const T& y1, const T& y2);

}


/*
 * Function definition.
 */

int orcs::utils::compare(double first, double second, double threshold) {
    return std::abs(first - second) < threshold ? 0 : (first < second ? -1 : 1);
}


template <class... T>
std::string orcs::utils::format(const std::string& str, T... args) {
    std::size_t buffer_size = std::max(static_cast<std::size_t>(1024), static_cast<std::size_t>(2 * str.length()));
    char* buffer = new char[buffer_size];
    sprintf(buffer, str.c_str(), args...);
    return std::string(buffer);
}


template <class T>
inline bool orcs::utils::overlap(const T& x1, const T& x2, const T& y1, const T& y2) {
    return (x1 <= y2 && y1 <= x2);
}


template <class T>
std::optional< std::tuple<T, T> > orcs::utils::intersection(const T& x1, const T& x2, const T& y1, const T& y2) {
    if (x1 <= y2 && y1 <= x2) {
        return std::make_tuple(std::max(x1, y1), std::min(x2, y2));
    }
    return {};
}

#endif

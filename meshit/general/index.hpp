#ifndef FILE_TEMPLATE_HPP
#define FILE_TEMPLATE_HPP

/**************************************************************************/
/* File:   template.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <iostream>
#include <unordered_map>

namespace meshit {
/**
  INDEX is a typedef for (at least) 4-byte integer
 */
typedef int INDEX;

class INDEX_2
{
 public:
    INDEX_2() { }

    INDEX_2(INDEX ai1, INDEX ai2)
    {
        i[0] = ai1;
        i[1] = ai2;
    }

    INDEX_2(const INDEX_2& in2)
    {
        i[0] = in2.i[0];
        i[1] = in2.i[1];
    }

    int operator==(const INDEX_2& in2) const { return i[0] == in2.i[0] && i[1] == in2.i[1]; }

    INDEX_2 Sort()
    {
        if (i[0] > i[1]) {
            std::swap(i[0], i[1]);
        }
        return *this;
    }

    INDEX& I1() { return i[0]; }
    INDEX& I2() { return i[1]; }
    const INDEX& I1() const { return i[0]; }
    const INDEX& I2() const { return i[1]; }

    int& operator[](size_t j) { return i[j]; }
    const int& operator[](size_t j) const { return i[j]; }

    friend std::ostream& operator<<(std::ostream& s, const INDEX_2& i2);

 protected:
    INDEX i[2];
};

template<class T> using INDEX_2_map = std::unordered_map<INDEX_2, T>;

}  // namespace meshit

namespace std {

template<>
struct hash<::meshit::INDEX_2>
{
    size_t operator()(const ::meshit::INDEX_2& idx) const { return idx.I1() ^ (idx.I2() << 16); }
};

}  // namespace std

#endif

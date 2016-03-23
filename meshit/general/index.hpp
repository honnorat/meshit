#ifndef FILE_TEMPLATE_HPP
#define FILE_TEMPLATE_HPP

/**************************************************************************/
/* File:   template.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <iostream>

namespace meshit
{
    /**
      INDEX is a typedef for (at least) 4-byte integer
     */
    typedef int INDEX;

    class INDEX_2;

    std::ostream& operator<<(std::ostream& s, const INDEX_2& i2);

    class INDEX_2
    {
        INDEX i[2];

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

        int operator==(const INDEX_2& in2) const
        {
            return i[0] == in2.i[0] && i[1] == in2.i[1];
        }

        INDEX_2 Sort()
        {
            if (i[0] > i[1]) {
                std::swap(i[0], i[1]);
            }
            return *this;
        }

        static INDEX_2 Sort(INDEX i1, INDEX i2)
        {
            if (i1 > i2)
                return INDEX_2(i2, i1);
            else
                return INDEX_2(i1, i2);
        }

        INDEX& I1()
        {
            return i[0];
        }

        INDEX& I2()
        {
            return i[1];
        }

        INDEX& I(size_t j)
        {
            return i[j - 1];
        }

        const INDEX& I1() const
        {
            return i[0];
        }

        const INDEX& I2() const
        {
            return i[1];
        }

        const INDEX& I(size_t j) const
        {
            return i[j - 1];
        }

        int& operator[](size_t j)
        {
            return i[j];
        }

        const int& operator[](size_t j) const
        {
            return i[j];
        }

        friend std::ostream& operator<<(std::ostream& s, const INDEX_2& i2);
    };

    inline INDEX_2 Sort(const INDEX_2& i2)
    {
        INDEX_2 tmp = i2;
        tmp.Sort();
        return tmp;
    }

    class INDEX_3
    {
        INDEX i[3];

     public:
        INDEX_3() { }

        INDEX_3(INDEX ai1, INDEX ai2, INDEX ai3)
        {
            i[0] = ai1;
            i[1] = ai2;
            i[2] = ai3;
        }

        INDEX_3(const INDEX_3& in2)
        {
            i[0] = in2.i[0];
            i[1] = in2.i[1];
            i[2] = in2.i[2];
        }

        static INDEX_3 Sort(INDEX_3 i3)
        {
            return i3.Sort();
        }

        static INDEX_3 Sort(int i1, int i2, int i3)
        {
            if (i1 > i2) std::swap(i1, i2);
            if (i2 > i3) std::swap(i2, i3);
            if (i1 > i2) std::swap(i1, i2);
            return INDEX_3(i1, i2, i3);
        }

        INDEX_3 Sort()
        {
            if (i[0] > i[1]) std::swap(i[0], i[1]);
            if (i[1] > i[2]) std::swap(i[1], i[2]);
            if (i[0] > i[1]) std::swap(i[0], i[1]);
            return *this;
        }

        int operator==(const INDEX_3& in2) const
        {
            return i[0] == in2.i[0] && i[1] == in2.i[1] && i[2] == in2.i[2];
        }

        INDEX& I1()
        {
            return i[0];
        }

        INDEX& I2()
        {
            return i[1];
        }

        INDEX& I3()
        {
            return i[2];
        }

        INDEX& I(size_t j)
        {
            return i[j - 1];
        }

        const INDEX& I1() const
        {
            return i[0];
        }

        const INDEX& I2() const
        {
            return i[1];
        }

        const INDEX& I3() const
        {
            return i[2];
        }

        const INDEX& I(size_t j) const
        {
            return i[j - 1];
        }

        int& operator[](size_t j)
        {
            return i[j];
        }

        const int& operator[](size_t j) const
        {
            return i[j];
        }

        friend std::ostream& operator<<(std::ostream& s, const INDEX_3& i3);
    };

    class INDEX_4
    {
        INDEX i[4];

     public:
        INDEX_4() { }

        INDEX_4(INDEX ai1, INDEX ai2, INDEX ai3, INDEX ai4)
        {
            i[0] = ai1;
            i[1] = ai2;
            i[2] = ai3;
            i[3] = ai4;
        }

        INDEX_4(const INDEX_4& in2)
        {
            i[0] = in2.i[0];
            i[1] = in2.i[1];
            i[2] = in2.i[2];
            i[3] = in2.i[3];
        }

        void Sort();

        int operator==(const INDEX_4& in2) const
        {
            return i[0] == in2.i[0] && i[1] == in2.i[1] &&
                   i[2] == in2.i[2] && i[3] == in2.i[3];
        }

        INDEX& I1()
        {
            return i[0];
        }

        INDEX& I2()
        {
            return i[1];
        }

        INDEX& I3()
        {
            return i[2];
        }

        INDEX& I4()
        {
            return i[3];
        }

        INDEX& I(int j)
        {
            return i[j - 1];
        }

        const INDEX& I1() const
        {
            return i[0];
        }

        const INDEX& I2() const
        {
            return i[1];
        }

        const INDEX& I3() const
        {
            return i[2];
        }

        const INDEX& I4() const
        {
            return i[3];
        }

        const INDEX& I(size_t j) const
        {
            return i[j - 1];
        }

        int& operator[](size_t j)
        {
            return i[j];
        }

        const int& operator[](size_t j) const
        {
            return i[j];
        }

        friend std::ostream& operator<<(std::ostream& s, const INDEX_4& i4);
    };

    inline bool operator<(const INDEX_4& a, const INDEX_4& b)
    {
        for (size_t j = 0; j < 4; j++) {
            if (a[j] < b[j]) return true;
            if (a[j] > b[j]) return false;
        }
        return false;
    }

}  // namespace meshit

#endif

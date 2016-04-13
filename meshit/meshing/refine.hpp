#ifndef REFINE_HPP
#define REFINE_HPP

#include "improve2.hpp"

namespace meshit
{
    class Refinement
    {
     public:
        Refinement() { }

        ~Refinement() { }

        void Refine(Mesh& mesh);

    };

}  // namespace meshit

#endif

#ifndef FEMTO_CONVENTION_HPP
#define FEMTO_CONVENTION_HPP

namespace femto {

// Nuclear convention: k cot(delta) = -1/a0 + re k^2/2
// Femto convention:  k cot(delta) =  1/f0 + d0 k^2/2  with f0 = -a0, d0 = re

inline double nuclearToFemtoF0(double a0) { return -a0; }
inline double femtoToNuclearA0(double f0) { return -f0; }
inline double nuclearToFemtoD0(double re) { return re; }

} // namespace femto

#endif

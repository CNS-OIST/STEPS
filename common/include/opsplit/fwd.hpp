#pragma once

#include <Omega_h_defines.hpp>

namespace zee {

namespace osh = Omega_h;

// Forward declaration
class Compdef;
class Diffdef;
class DistMesh;
class SReacdef;
class Reacdef;
class Patchdef;

template <typename RNG>
class Simulation;
class Statedef;
struct ContainerKProcSegment;
class KProcState;

enum class KProcType : unsigned {
    Reac = 0x10,
    Diff = 0x8,
    SReac = 0x4,
    VDepSReac = 0x2,
    SDiff = 0x1
};

inline bool any(KProcType type) {
    return static_cast<unsigned>(type) != 0;
}

constexpr inline KProcType operator|(KProcType lhs, KProcType rhs) {
    return static_cast<KProcType>(static_cast<unsigned>(lhs) | static_cast<unsigned>(rhs));
}

constexpr inline KProcType operator&(KProcType lhs, KProcType rhs) {
    return static_cast<KProcType>(static_cast<unsigned>(lhs) & static_cast<unsigned>(rhs));
}

constexpr inline KProcType operator~(KProcType rhs) {
    return static_cast<KProcType>(~static_cast<unsigned>(rhs));
}

static const auto ALL_KPROC_TYPES = KProcType::Reac | KProcType::Diff | KProcType::SReac |
                                    KProcType::VDepSReac | KProcType::SDiff;

}  // namespace zee

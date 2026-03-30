#pragma once

namespace steps::dist::kproc {

class KProcState;

enum class KProcType : unsigned {
    Reac = 0,
    Diff = 1,
    SReac = 2,
    VDepSReac = 3,
    GHKSReac = 4,
    ComplexReac = 5,
    ComplexSReac = 6,
    VDepComplexSReac = 7,
    ComplexGHKSReac = 8,
};
constexpr unsigned num_kproc_types() {
    return 9;
}

}  // namespace steps::dist::kproc

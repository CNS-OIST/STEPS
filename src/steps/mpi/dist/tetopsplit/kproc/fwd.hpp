#pragma once

namespace steps::dist::kproc {

class KProcState;

enum class KProcType : unsigned { Reac = 0, Diff = 1, SReac = 2, VDepSReac = 3, GHKSReac = 4 };
constexpr unsigned num_kproc_types() {
    return 5;
}

}  // namespace steps::dist::kproc

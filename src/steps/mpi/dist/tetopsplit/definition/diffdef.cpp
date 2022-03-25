#include "diffdef.hpp"

#include "compdef.hpp"
#include "statedef.hpp"


namespace steps {
namespace dist {

Diffdef::Diffdef(const Compdef &compdef, container::kproc_id kproc,
                 container::diffusion_id t_diffusion,
                 container::species_id species, osh::Real t_dcst)
    : pCompdef(compdef), kprocContainerIdx(kproc), diffusion(t_diffusion),
      specContainerIdx(species), dcst(t_dcst) {}

void Diffdef::report(std::ostream& ostr) const {
    ostr << "Diffusion Report" << std::endl;
    ostr << "KProc Container Idx: " << kprocContainerIdx << " Diff Container Idx: " << diffusion
         << std::endl;
    const auto spec_model_idx = pCompdef.getSpecModelIdx(specContainerIdx);
    ostr << "Diffusion Spec :" << pCompdef.statedef().getSpecID(spec_model_idx)
         << " Model Idx: " << spec_model_idx << " Container Idx: " << specContainerIdx
         << " DCST: " << dcst << std::endl;
}

}  // namespace dist
}  // namespace steps
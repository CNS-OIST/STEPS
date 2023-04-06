#pragma once

#include "geom/dist/distmesh.hpp"
#include "solver/api.hpp"
#include "definition/statedef.hpp"
#include "simulation.hpp"
#include "fwd.hpp"

namespace steps {
namespace dist {

class TetOpSplitBase {
  public:
    TetOpSplitBase(){};
    virtual ~TetOpSplitBase(){};

    virtual std::string getSolverName() const = 0;
    virtual std::string getSolverDesc() const = 0;
    virtual std::string getSolverAuthors() const = 0;
    virtual std::string getSolverEmail() const = 0;

    virtual void reset() = 0;
    virtual void run(double seconds) = 0;
    virtual double getTime() const = 0;

    virtual double getCompCount(std::string const &c,
                                std::string const &s) const = 0;
    virtual void setCompCount(std::string const &c, std::string const &s,
                              double n,
                              const math::DistributionMethod distribution) = 0;
    virtual double getCompConc(std::string const &c,
                               std::string const &s) const = 0;
    virtual void setCompConc(std::string const &c, std::string const &s,
                             double conc,
                             const math::DistributionMethod distribution) = 0;
    virtual double getPatchCount(std::string const &p,
                                 std::string const &s) const = 0;
    virtual void setPatchCount(std::string const &p, std::string const &s,
                               double n,
                               const math::DistributionMethod distribution) = 0;
    virtual void setPatchSReacK(std::string const &p, std::string const &r,
                                double kf) = 0;

    virtual double getPatchMaxV(std::string const &p) const = 0;

    virtual void setMembPotential(const std::string &memb, double pot) = 0;
    virtual double getVertIClamp(osh::GO vert, bool local) const = 0;
    virtual void setVertIClamp(osh::GO vert, double current, bool local) = 0;

    virtual double getTetCount(osh::GO tet, const std::string &s,
                               bool local) const = 0;
    virtual double getTetConc(osh::GO tet, const std::string &s,
                              bool local) const = 0;
    virtual void setTetCount(osh::GO tet, const std::string &s,
                             double count, bool local) = 0;
    virtual void setTetConc(osh::GO tet, const std::string &s, double conc, bool local) = 0;
    virtual double getTriCount(osh::GO tri, const std::string &s,
                               bool local) const = 0;
    virtual void setTriCount(osh::GO tri, const std::string &s,
                             double count, bool local) = 0;
    virtual double getVertV(osh::GO vert, bool local) const = 0;
    virtual double getTriV(osh::GO tet, bool local) const = 0;
    virtual double getTetV(osh::GO tet, bool local) const = 0;
    virtual double getTriOhmicI(osh::GO tet, const std::string &curr,
                                bool local) const = 0;
    virtual double getTriGHKI(osh::GO tet, const std::string &curr,
                              bool local) const = 0;

    virtual std::vector<double>
    getBatchTetCounts(const std::vector<osh::GO> &tets, std::string const &s,
                      bool local) const = 0;
    virtual std::vector<double>
    getBatchTetConcs(const std::vector<osh::GO> &tets, std::string const &s,
                     bool local) const = 0;
    virtual void setBatchTetCounts(const std::vector<osh::GO> &tets,
                                   std::string const &s,
                                   std::vector<double> &counts, bool local) = 0;
    virtual void setBatchTetConcs(const std::vector<osh::GO> &tets,
                                  std::string const &s,
                                  std::vector<double> &concs, bool local) = 0;
    virtual std::vector<double>
    getBatchTriCounts(const std::vector<osh::GO> &tris, std::string const &s,
                      bool local) const = 0;
    virtual void setBatchTriCounts(const std::vector<osh::GO> &tris,
                                   std::string const &s,
                                   std::vector<double> &counts, bool local) = 0;
    virtual void getBatchTetCountsNP(const osh::GO *indices, size_t input_size,
                                     std::string const &s, double *counts,
                                     size_t output_size, bool local) const = 0;
    virtual void getBatchTetConcsNP(const osh::GO *indices, size_t input_size,
                                    std::string const &s, double *concs,
                                    size_t output_size, bool local) const = 0;
    virtual void setBatchTetCountsNP(const osh::GO *indices, size_t input_size,
                                     std::string const &s, double *counts,
                                     size_t output_size, bool local) = 0;
    virtual void setBatchTetConcsNP(const osh::GO *indices, size_t input_size,
                                    std::string const &s, double *concs,
                                    size_t output_size, bool local) = 0;
    virtual void getBatchTriCountsNP(const osh::GO *indices, size_t input_size,
                                     std::string const &s, double *counts,
                                     size_t output_size, bool local) const = 0;
    virtual void setBatchTriCountsNP(const osh::GO *indices, size_t input_size,
                                     std::string const &s, double *counts,
                                     size_t output_size, bool local) = 0;
    virtual void getBatchVertVsNP(const osh::GO *indices, size_t input_size,
                                  double *voltages, size_t output_size,
                                  bool local) const = 0;
    virtual void getBatchTriVsNP(const osh::GO *indices, size_t input_size,
                                 double *voltages, size_t output_size,
                                 bool local) const = 0;
    virtual void getBatchTetVsNP(const osh::GO *indices, size_t input_size,
                                 double *voltages, size_t output_size,
                                 bool local) const = 0;
    virtual void getBatchTriOhmicIsNP(const osh::GO *indices, size_t input_size,
                                      const std::string &curr, double *currents,
                                      size_t output_size, bool local) const = 0;
    virtual void getBatchTriGHKIsNP(const osh::GO *indices, size_t input_size,
                                    const std::string &curr, double *currents,
                                    size_t output_size, bool local) const = 0;

    virtual void setDiffBoundaryDiffusionActive(const std::string &name,
                                                const std::string &spec,
                                                bool active) = 0;
    virtual bool getDiffBoundaryDiffusionActive(const std::string &name,
                                                const std::string &spec) = 0;
    virtual void setDiffApplyThreshold(int threshold) = 0;
    virtual void setMembIClamp(const std::string &memb, double stim) = 0;

    /// Get temperature
    virtual osh::Real getTemp() const noexcept = 0;
    /// Set temperature
    virtual void setTemp(const osh::Real temp) noexcept = 0;
#if USE_PETSC
    virtual double getEfieldDT() const = 0;
    virtual void setEfieldDT(double dt) = 0;
    virtual void setEfieldTolerances(double atol, double rtol,
                                     KSPNormType norm_type) = 0;
#endif // USE_PETSC
    virtual double getCompTime() const noexcept = 0;
    virtual double getSyncTime() const noexcept = 0;
    virtual double getIdleTime() const noexcept = 0;
    virtual double getEFieldTime() const noexcept = 0;
    virtual double getRDTime() const noexcept = 0;
    virtual double getDiffusionTime() const noexcept = 0;
    virtual double getReactionTime() const noexcept = 0;
    virtual double getDataExchangeTime() const noexcept = 0;
    virtual unsigned long long getDiffExtent(bool local) const = 0;
    virtual unsigned long long getReacExtent(bool local) const = 0;
};

template <steps::dist::SSAMethod SSA = steps::dist::SSAMethod::SSA,
          steps::dist::NextEventSearchMethod SearchMethod =
              steps::dist::NextEventSearchMethod::Direct>
class TetOpSplit : public TetOpSplitBase {
  public:
    TetOpSplit(steps::model::Model &model, steps::dist::DistMesh &mesh,
               const rng::RNGptr &r, bool indepKProcs);
    virtual ~TetOpSplit();

    /**
     * \name Solver information
     * \{
     */

    std::string getSolverName() const override;
    std::string getSolverDesc() const override;
    std::string getSolverAuthors() const override;
    std::string getSolverEmail() const override;

    /**
     * \}
     */

    /**
     * \name Solver general controls
     * \{
     */

    void reset() override;
    void run(double seconds) override;
    double getTime() const override;

    /**
     * \}
     */

    /**
     * \name Solver compartment controls
     * \{
     */

    double getCompCount(std::string const &c,
                        std::string const &s) const override;
    void setCompCount(std::string const &c, std::string const &s, double n,
                      const math::DistributionMethod distribution) override;

    double getCompConc(std::string const &c,
                       std::string const &s) const override;
    void setCompConc(std::string const &c, std::string const &s, double conc,
                     const math::DistributionMethod distribution) override;

    /**
     * \}
     */

    /**
     * \name Solver patch controls
     * \{
     */

    double getPatchCount(std::string const &p,
                         std::string const &s) const override;

    void setPatchCount(std::string const &p, std::string const &s, double n,
                       const math::DistributionMethod distribution) override;

    void setPatchSReacK(std::string const &p, std::string const &r,
                        double kf) override;

    double getPatchMaxV(std::string const &p) const override;

    /**
     * \}
     */

    /**
     * \name Solver membrane controls
     * \{
     */

    void setMembPotential(const std::string &memb, double pot) override;

    /**
     * \}
     */

    /**
     * \name Single element data access
     * \{
     */

    double getVertIClamp(const osh::GO vert, bool local) const override;

    void setVertIClamp(const osh::GO vert, const double value,
                       bool local) override;

    double getTetCount(osh::GO tet, const std::string &s,
                       bool local) const override;

    double getTetConc(osh::GO tet, const std::string &s,
                      bool local) const override;

    void setTetCount(osh::GO tet, const std::string &s, double count,
                     bool local) override;

    void setTetConc(osh::GO tet, const std::string &s, double conc,
                    bool local) override;

    double getTriCount(osh::GO tri, const std::string &s,
                       bool local) const override;

    void setTriCount(osh::GO tri, const std::string &s, double count,
                     bool local) override;

    double getVertV(osh::GO vert, bool local) const override;

    double getTriV(osh::GO tri, bool local) const override;

    double getTetV(osh::GO tet, bool local) const override;

    double getTriOhmicI(osh::GO tet, const std::string &curr,
                        bool local) const override;

    double getTriGHKI(osh::GO tet, const std::string &curr,
                      bool local) const override;

    /**
     * \}
     */

    /**
     * \name Batch data access
     * \{
     */

    std::vector<double> getBatchTetCounts(const std::vector<osh::GO> &tets,
                                          std::string const &s,
                                          bool local) const override;

    std::vector<double> getBatchTetConcs(const std::vector<osh::GO> &tets,
                                         std::string const &s,
                                         bool local) const override;

    void setBatchTetCounts(const std::vector<osh::GO> &tets,
                           std::string const &s, std::vector<double> &counts,
                           bool local) override;

    void setBatchTetConcs(const std::vector<osh::GO> &tets,
                          std::string const &s, std::vector<double> &concs,
                          bool local) override;

    std::vector<double> getBatchTriCounts(const std::vector<osh::GO> &tris,
                                          std::string const &s,
                                          bool local) const override;

    void setBatchTriCounts(const std::vector<osh::GO> &tris,
                           std::string const &s, std::vector<double> &counts,
                           bool local) override;

    void getBatchTetCountsNP(const osh::GO *indices, size_t input_size,
                             std::string const &s, double *counts,
                             size_t output_siz, bool local) const override;

    void getBatchTetConcsNP(const osh::GO *indices, size_t input_size,
                            std::string const &s, double *concs,
                            size_t output_size, bool local) const override;

    void setBatchTetCountsNP(const osh::GO *indices, size_t input_size,
                             std::string const &s, double *counts,
                             size_t output_size, bool local) override;

    void setBatchTetConcsNP(const osh::GO *indices, size_t input_size,
                            std::string const &s, double *concs,
                            size_t output_size, bool local) override;

    void getBatchTriCountsNP(const osh::GO *indices, size_t input_size,
                             std::string const &s, double *counts,
                             size_t output_size, bool local) const override;

    void setBatchTriCountsNP(const osh::GO *indices, size_t input_size,
                             std::string const &s, double *counts,
                             size_t output_size, bool local) override;

    void getBatchVertVsNP(const osh::GO *indices, size_t input_size,
                          double *voltages, size_t output_size,
                          bool local) const override;

    void getBatchTriVsNP(const osh::GO *indices, size_t input_size,
                         double *voltages, size_t output_size,
                         bool local) const override;

    void getBatchTetVsNP(const osh::GO *indices, size_t input_size,
                         double *voltages, size_t output_size,
                         bool local) const override;

    void getBatchTriOhmicIsNP(const osh::GO *indices, size_t input_size,
                              const std::string &curr, double *currents,
                              size_t output_size, bool local) const override;

    void getBatchTriGHKIsNP(const osh::GO *indices, size_t input_size,
                            const std::string &curr, double *currents,
                            size_t output_size, bool local) const override;

    /**
     * \}
     */

    void setDiffBoundaryDiffusionActive(const std::string &name,
                                        const std::string &spec,
                                        bool active) override;

    bool getDiffBoundaryDiffusionActive(const std::string &name,
                                        const std::string &spec) override;

    /**
     * \name Counters and thresholds
     * \{
     */

    void setDiffApplyThreshold(int threshold) override;

    /**
     * \}
     */
    void setMembIClamp(const std::string &memb, double stim) override;

    /// Get temperature
    osh::Real getTemp() const noexcept override { return sim->getTemp(); }
    /// Set temperature
    void setTemp(const osh::Real temp) noexcept override { sim->setTemp(temp); }

#if USE_PETSC
    double getEfieldDT() const override;
    void setEfieldDT(double dt) override;
    void setEfieldTolerances(double atol, double rtol,
                             KSPNormType norm_type) override;
#endif // USE_PETSC

    /**
     * \name Monitoring methods
     * \{
     */

    double getCompTime() const noexcept override { return comptime; }
    double getSyncTime() const noexcept override { return synctime; }
    double getIdleTime() const noexcept override { return idletime; }
    double getEFieldTime() const noexcept override {
        return sim->getElapsedEField();
    }
    double getRDTime() const noexcept override {
        return sim->getElapsedDiff() + sim->getElapsedSSA();
    }
    double getDiffusionTime() const noexcept override {
        return sim->getElapsedDiff();
    }
    double getReactionTime() const noexcept override {
        return sim->getElapsedSSA();
    }
    double getDataExchangeTime() const noexcept override {
        return dataexchangetime;
    }

    unsigned long long getDiffExtent(bool local) const override {
        return sim->getDiffOpExtent(local);
    }
    unsigned long long getReacExtent(bool local) const override {
        return sim->getSSAOpExtent(local);
    }

    /**
     * \}
     */

  private:
    steps::dist::DistMesh& meshref;
    std::unique_ptr<steps::dist::OmegaHSimulation<SSA, steps::rng::RNG,
                                                  osh::I64, SearchMethod>>
        sim;

    double comptime{};
    double synctime{};
    double idletime{};
    double efieldtime{};
    double rdtime{};
    double dataexchangetime{};
};

}  // namespace dist
}  // namespace steps

#pragma once

#include "definition/statedef.hpp"
#include "fwd.hpp"
#include "geom/dist/distmesh.hpp"
#include "simulation.hpp"
#include "solver/api.hpp"

namespace steps::dist {

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

    virtual double getCompSpecCount(std::string const& c, std::string const& s) const = 0;
    virtual void setCompSpecCount(std::string const& c,
                                  std::string const& s,
                                  double n,
                                  const math::DistributionMethod distribution) = 0;
    virtual double getCompSpecConc(std::string const& c, std::string const& s) const = 0;
    virtual void setCompSpecConc(std::string const& c,
                                 std::string const& s,
                                 double conc,
                                 const math::DistributionMethod distribution) = 0;
    virtual double getPatchSpecCount(std::string const& p, std::string const& s) const = 0;
    virtual void setPatchSpecCount(std::string const& p,
                                   std::string const& s,
                                   double n,
                                   const math::DistributionMethod distribution) = 0;
    virtual void setPatchSReacK(std::string const& p, std::string const& r, double kf) = 0;

    virtual double getPatchMaxV(std::string const& p) const = 0;

    virtual void setMembPotential(const std::string& memb, double pot) = 0;
    virtual void setMembRes(const std::string& membrane,
                            osh::Real resistivity,
                            osh::Real reversal_potential) = 0;
    virtual MembraneResistivity getMembRes(const std::string& membrane) const = 0;
    virtual MembraneResistivity getTriRes(osh::GO tri, bool local) const = 0;
    virtual void setTriRes(const osh::GO tri, osh::Real res, osh::Real erev, bool local) = 0;
    virtual double getTriCapac(osh::GO tri, bool local) const = 0;
    virtual void setTriCapac(const osh::GO tri, osh::Real c, bool local) = 0;
    virtual double getVertIClamp(osh::GO vert, bool local) const = 0;
    virtual void setVertIClamp(osh::GO vert, double current, bool local) = 0;

#if USE_PETSC
    virtual double getTriOhmicErev(osh::GO triangle,
                                   const std::string& ohmic_current,
                                   bool local) const = 0;
    virtual void setTriOhmicErev(osh::GO triangle,
                                 const std::string& ohmic_current,
                                 double reversal_potential,
                                 bool local) = 0;
    virtual void getBatchTriOhmicErevsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const std::string& ohmic_current,
                                         double* rv,
                                         size_t output_size,
                                         bool local) const = 0;
#endif  // USE_PETSC

    virtual double getTetSpecCount(osh::GO tet, const std::string& s, bool local) const = 0;
    virtual double getTetSpecConc(osh::GO tet, const std::string& s, bool local) const = 0;
    virtual void setTetSpecCount(osh::GO tet, const std::string& s, double count, bool local) = 0;
    virtual void setTetSpecConc(osh::GO tet, const std::string& s, double conc, bool local) = 0;
    virtual double getTriSpecCount(osh::GO tri, const std::string& s, bool local) const = 0;
    virtual void setTriSpecCount(osh::GO tri, const std::string& s, double count, bool local) = 0;
    virtual double getVertV(osh::GO vert, bool local) const = 0;
    virtual double getTriV(osh::GO tet, bool local) const = 0;
    virtual double getTetV(osh::GO tet, bool local) const = 0;
    virtual double getTriOhmicI(osh::GO tet, const std::string& curr, bool local) const = 0;
    virtual double getTriGHKI(osh::GO tet, const std::string& curr, bool local) const = 0;

    virtual std::vector<double> getBatchTetSpecCounts(const std::vector<osh::GO>& tets,
                                                      std::string const& s,
                                                      bool local) const = 0;
    virtual std::vector<double> getBatchTetSpecConcs(const std::vector<osh::GO>& tets,
                                                     std::string const& s,
                                                     bool local) const = 0;
    virtual void setBatchTetSpecCounts(const std::vector<osh::GO>& tets,
                                       std::string const& s,
                                       std::vector<double>& counts,
                                       bool local) = 0;
    virtual void setBatchTetSpecConcs(const std::vector<osh::GO>& tets,
                                      std::string const& s,
                                      std::vector<double>& concs,
                                      bool local) = 0;
    virtual std::vector<double> getBatchTriSpecCounts(const std::vector<osh::GO>& tris,
                                                      std::string const& s,
                                                      bool local) const = 0;
    virtual void setBatchTriSpecCounts(const std::vector<osh::GO>& tris,
                                       std::string const& s,
                                       std::vector<double>& counts,
                                       bool local) = 0;
    virtual void getBatchTetSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         std::string const& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) const = 0;
    virtual void getBatchTetSpecConcsNP(const osh::GO* indices,
                                        size_t input_size,
                                        std::string const& s,
                                        double* concs,
                                        size_t output_size,
                                        bool local) const = 0;
    virtual void setBatchTetSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         std::string const& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) = 0;
    virtual void setBatchTetSpecConcsNP(const osh::GO* indices,
                                        size_t input_size,
                                        std::string const& s,
                                        double* concs,
                                        size_t output_size,
                                        bool local) = 0;
    virtual void getBatchTriSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         std::string const& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) const = 0;
    virtual void setBatchTriSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         std::string const& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) = 0;
    virtual void getBatchVertVsNP(const osh::GO* indices,
                                  size_t input_size,
                                  double* voltages,
                                  size_t output_size,
                                  bool local) const = 0;
    virtual void getBatchTriVsNP(const osh::GO* indices,
                                 size_t input_size,
                                 double* voltages,
                                 size_t output_size,
                                 bool local) const = 0;
    virtual void getBatchTetVsNP(const osh::GO* indices,
                                 size_t input_size,
                                 double* voltages,
                                 size_t output_size,
                                 bool local) const = 0;
    virtual void getBatchTriOhmicIsNP(const osh::GO* indices,
                                      size_t input_size,
                                      const std::string& curr,
                                      double* currents,
                                      size_t output_size,
                                      bool local) const = 0;
    virtual void getBatchTriGHKIsNP(const osh::GO* indices,
                                    size_t input_size,
                                    const std::string& curr,
                                    double* currents,
                                    size_t output_size,
                                    bool local) const = 0;

    virtual void setDiffBoundarySpecDiffusionActive(const std::string& name,
                                                    const std::string& spec,
                                                    bool active) = 0;
    virtual bool getDiffBoundarySpecDiffusionActive(const std::string& name,
                                                    const std::string& spec) = 0;
    virtual void setDiffApplyThreshold(int threshold) = 0;
    virtual void setMembIClamp(const std::string& memb, double stim) = 0;

    /// Get temperature
    virtual osh::Real getTemp() const noexcept = 0;
    /// Set temperature
    virtual void setTemp(const osh::Real temp) noexcept = 0;
#if USE_PETSC
    virtual double getEfieldDT() const = 0;
    virtual void setEfieldDT(double dt) = 0;
    virtual void setPetscOptions(const std::string& s) = 0;
#endif  // USE_PETSC
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
    virtual void dumpDepGraphToFile(const std::string& path) const = 0;
};

template <SSAMethod SSA = SSAMethod::SSA,
          NextEventSearchMethod SearchMethod = NextEventSearchMethod::Direct>
class TetOpSplit: public TetOpSplitBase {
  public:
    TetOpSplit(steps::model::Model& model,
               steps::dist::DistMesh& mesh,
               const rng::RNGptr& r,
               bool indepKProcs,
               bool isEfield);
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

    double getCompSpecCount(std::string const& c, std::string const& s) const override;
    void setCompSpecCount(std::string const& c,
                          std::string const& s,
                          double n,
                          const math::DistributionMethod distribution) override;

    double getCompSpecConc(std::string const& c, std::string const& s) const override;
    void setCompSpecConc(std::string const& c,
                         std::string const& s,
                         double conc,
                         const math::DistributionMethod distribution) override;

    /**
     * \}
     */

    /**
     * \name Solver patch controls
     * \{
     */

    double getPatchSpecCount(std::string const& p, std::string const& s) const override;

    void setPatchSpecCount(std::string const& p,
                           std::string const& s,
                           double n,
                           const math::DistributionMethod distribution) override;

    void setPatchSReacK(std::string const& p, std::string const& r, double kf) override;

    double getPatchMaxV(std::string const& p) const override;

    /**
     * \}
     */

    /**
     * \name Solver membrane controls
     * \{
     */

    void setMembPotential(const std::string& memb, double pot) override;

    void setMembRes(const std::string& membrane,
                    osh::Real resistivity,
                    osh::Real reversal_potential) override;

    MembraneResistivity getMembRes(const std::string& membrane) const override;

    MembraneResistivity getTriRes(osh::GO tri, bool local) const override;

    void setTriRes(const osh::GO tri, osh::Real res, osh::Real erev, bool local) override;

    /**
     * \}
     */

    /**
     * \name Single element data access
     * \{
     */

    double getTriCapac(osh::GO tri, bool local) const override;

    void setTriCapac(const osh::GO tri, osh::Real c, bool local) override;

    double getVertIClamp(const osh::GO vert, bool local) const override;

    void setVertIClamp(const osh::GO vert, const double value, bool local) override;

#if USE_PETSC
    double getTriOhmicErev(osh::GO triangle,
                           const std::string& ohmic_current,
                           bool local) const override;
    void setTriOhmicErev(osh::GO triangle,
                         const std::string& ohmic_current,
                         double reversal_potential,
                         bool local) override;
    void getBatchTriOhmicErevsNP(const osh::GO* indices,
                                 size_t input_size,
                                 const std::string& ohmic_current,
                                 double* rv,
                                 size_t output_size,
                                 bool local) const override;

#endif  // USE_PETSC

    double getTetSpecCount(osh::GO tet, const std::string& s, bool local) const override;

    double getTetSpecConc(osh::GO tet, const std::string& s, bool local) const override;

    void setTetSpecCount(osh::GO tet, const std::string& s, double count, bool local) override;

    void setTetSpecConc(osh::GO tet, const std::string& s, double conc, bool local) override;

    double getTriSpecCount(osh::GO tri, const std::string& s, bool local) const override;

    void setTriSpecCount(osh::GO tri, const std::string& s, double count, bool local) override;

    double getVertV(osh::GO vert, bool local) const override;

    double getTriV(osh::GO tri, bool local) const override;

    double getTetV(osh::GO tet, bool local) const override;

    double getTriOhmicI(osh::GO tet, const std::string& curr, bool local) const override;

    double getTriGHKI(osh::GO tet, const std::string& curr, bool local) const override;

    /**
     * \}
     */

    /**
     * \name Batch data access
     * \{
     */

    std::vector<double> getBatchTetSpecCounts(const std::vector<osh::GO>& tets,
                                              std::string const& s,
                                              bool local) const override;

    std::vector<double> getBatchTetSpecConcs(const std::vector<osh::GO>& tets,
                                             std::string const& s,
                                             bool local) const override;

    void setBatchTetSpecCounts(const std::vector<osh::GO>& tets,
                               std::string const& s,
                               std::vector<double>& counts,
                               bool local) override;

    void setBatchTetSpecConcs(const std::vector<osh::GO>& tets,
                              std::string const& s,
                              std::vector<double>& concs,
                              bool local) override;

    std::vector<double> getBatchTriSpecCounts(const std::vector<osh::GO>& tris,
                                              std::string const& s,
                                              bool local) const override;

    void setBatchTriSpecCounts(const std::vector<osh::GO>& tris,
                               std::string const& s,
                               std::vector<double>& counts,
                               bool local) override;

    void getBatchTetSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 std::string const& s,
                                 double* counts,
                                 size_t output_siz,
                                 bool local) const override;

    void getBatchTetSpecConcsNP(const osh::GO* indices,
                                size_t input_size,
                                std::string const& s,
                                double* concs,
                                size_t output_size,
                                bool local) const override;

    void setBatchTetSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 std::string const& s,
                                 double* counts,
                                 size_t output_size,
                                 bool local) override;

    void setBatchTetSpecConcsNP(const osh::GO* indices,
                                size_t input_size,
                                std::string const& s,
                                double* concs,
                                size_t output_size,
                                bool local) override;

    void getBatchTriSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 std::string const& s,
                                 double* counts,
                                 size_t output_size,
                                 bool local) const override;

    void setBatchTriSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 std::string const& s,
                                 double* counts,
                                 size_t output_size,
                                 bool local) override;

    void getBatchVertVsNP(const osh::GO* indices,
                          size_t input_size,
                          double* voltages,
                          size_t output_size,
                          bool local) const override;

    void getBatchTriVsNP(const osh::GO* indices,
                         size_t input_size,
                         double* voltages,
                         size_t output_size,
                         bool local) const override;

    void getBatchTetVsNP(const osh::GO* indices,
                         size_t input_size,
                         double* voltages,
                         size_t output_size,
                         bool local) const override;

    void getBatchTriOhmicIsNP(const osh::GO* indices,
                              size_t input_size,
                              const std::string& curr,
                              double* currents,
                              size_t output_size,
                              bool local) const override;

    void getBatchTriGHKIsNP(const osh::GO* indices,
                            size_t input_size,
                            const std::string& curr,
                            double* currents,
                            size_t output_size,
                            bool local) const override;
    /**
     * \}
     */

    void setDiffBoundarySpecDiffusionActive(const std::string& name,
                                            const std::string& spec,
                                            bool active) override;

    bool getDiffBoundarySpecDiffusionActive(const std::string& name,
                                            const std::string& spec) override;

    /**
     * \name Counters and thresholds
     * \{
     */

    void setDiffApplyThreshold(int threshold) override;

    /**
     * \}
     */
    void setMembIClamp(const std::string& memb, double stim) override;

    /// Get temperature
    osh::Real getTemp() const noexcept override {
        return sim->getTemp();
    }
    /// Set temperature
    void setTemp(const osh::Real temp) noexcept override {
        sim->setTemp(temp);
    }

#if USE_PETSC
    double getEfieldDT() const override;
    void setEfieldDT(double dt) override;
    void setPetscOptions(const std::string& s) override;
#endif  // USE_PETSC

    /**
     * \name Monitoring methods
     * \{
     */

    double getCompTime() const noexcept override {
        return comptime;
    }
    double getSyncTime() const noexcept override {
        return synctime;
    }
    double getIdleTime() const noexcept override {
        return idletime;
    }
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

    /// Dump the dependency graph of kproc in a file specified by path
    void dumpDepGraphToFile(const std::string& path) const override;

    /**
     * \}
     */

  private:
    steps::dist::DistMesh& meshref;
    std::unique_ptr<steps::dist::OmegaHSimulation<SSA, SearchMethod>> sim;

    double comptime{};
    double synctime{};
    double idletime{};
    double dataexchangetime{};
};

}  // namespace steps::dist

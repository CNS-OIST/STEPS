/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

#pragma once

#include <numeric>
#include <cstddef>
#include <algorithm>

namespace steps {
namespace math {

/** Reservoir-based order samping with adjusted Pareto orders.
 *
 * Sample n items from a given iterator range with inclusion probabilities
 * pi[i].
 * 
 * Parameter is constructed from n and supplied pi; N is taken as the
 * length of pi. Template parameter weight_type is expected to be a
 * floating point type.
 *
 * Adjusted Pareto sampling selects the lowest ranked items where the rank is
 * given by Q[i] = U[i]/(1-U[i]) . q[i], where the U[i] are uniformly
 * distributed on [0,1) and the q[i] are given by:
 *
 *     q[i] = (1-pi[i])/pi[i] . exp(pi[i](1-pi[i])(pi[i]-1/2)/d^2)
 *        d = sum p[i](1-p[i])
 *
 *  In this implementation, the exponential is approximated by a
 *  second order polynomial, as the argument is small.
 */

template <typename weight_type=double,typename size_type=std::size_t>
struct adjusted_pareto_sampler {
    adjusted_pareto_sampler() {}

    template <typename Iter>
    adjusted_pareto_sampler(size_t n,Iter pi_begin,Iter pi_end):
        P(n,pi_begin,pi_end) {}

    struct param_type {
        size_type n=0;
        std::vector<weight_type> qcoef;

        param_type() {}
        template <typename Iter>
        param_type(size_type n_,Iter pi_begin,Iter pi_end):
            n(n_), qcoef(pi_begin,pi_end)
        {
            using fp=weight_type;

            // Normalize supplied pi so that sum = n, and compute d
            
            fp pi_sum=std::accumulate(qcoef.begin(),qcoef.end(),fp(0));
            fp rescale=n/pi_sum;

            fp d=0;
            for (auto &p: qcoef) {
                p*=rescale;
                if (p>1) p=1;
                d+=p*(1-p);
            }

            // Coefficients q in order sampling (see above)

            fp ood2=1/(d*d);
            for (auto &q: qcoef) {
                fp loga=q*(1-q)*(q-fp(0.5))*ood2;
                fp a=1+loga+fp(0.5)*loga*loga;
                q=(1-q)/q*a;
            }
        }
    } P;

    // Minium population size to sample.
    size_type size() const { return P.n; }

    // (Inclusive) min and max sample sizes
    size_type min() const { return P.n; }
    size_type max() const { return P.n; }

    void param(const param_type &P_) { P=P_; }
    const param_type &param() const { return P; }

    template <typename InIter,typename OutRAIter,typename Rng>
    size_type operator()(InIter b,InIter e,OutRAIter o,Rng &g) {
        std::uniform_real_distribution<weight_type> U;

        using key=std::pair<weight_type,size_type>;
        std::vector<key> heap;
        heap.reserve(size());

        size_type i=0;
        for (; i<P.n && b!=e; ++i) {
            auto u=U(g);
            heap.emplace_back(P.qcoef[i]*u/(1-u),i);
            o[i]=*b++;
        }

        if (i<P.n) return i; // short population

        std::make_heap(heap.begin(),heap.end());
        std::size_t N=P.qcoef.size();
        for (; i<N && b!=e; ++b) {
            auto u=U(g);
            weight_type q=P.qcoef[i++]*u/(1-u);

            if (q<heap.front().first) {
                key k{q,heap.front().second};
                std::pop_heap(heap.begin(),heap.end());
                heap.back()=k;
                o[k.second]=*b;
                std::push_heap(heap.begin(),heap.end());
            }
        }

        return P.n;
    }
};


} // namespace math
} // namespace steps

/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

#include "memb.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_set>

#include <easylogging++.h>

#include "tmcomp.hpp"
#include "tmpatch.hpp"

#include "util/collections.hpp"
#include "util/error.hpp"

namespace steps {
namespace tetmesh {

Memb::Memb(std::string id, Tetmesh *container,
           std::vector<TmPatch *> const &patches, bool verify, uint opt_method,
           double search_percent, std::string opt_file_name)
    : pID(std::move(id)), pTetmesh(container), pOpt_method(opt_method),
      pOpt_file_name(std::move(opt_file_name)),
      pSearch_percent(search_percent) {
  using steps::ArgErr;

  ArgErrLogIf(pTetmesh == nullptr,
              "No mesh provided to Membrane initializer function.");
  ArgErrLogIf(patches.empty(),
              "No Patches provided to Membrane initializer function.");
  ArgErrLogIf(pOpt_method != 1 && pOpt_method != 2,
              "Unknown optimization method. Choices are 1 or 2.");
  ArgErrLogIf(pSearch_percent > 100.0,
              "Search percentage is greater than 100.");
  ArgErrLogIf(pSearch_percent <= 0.0,
              "Search percentage must be greater than 0.");

  std::set<tetrahedron_id_t> tets_set;
  for (auto p : patches) {
    const auto &patch_tris = p->_getAllTriIndices();
    pTri_indices.insert(pTri_indices.end(), patch_tris.begin(),
                        patch_tris.end());

    auto comp = dynamic_cast<TmComp *>(p->getIComp());
    ArgErrLogIf(comp == nullptr,
                "Conduction volume (inner compartment(s) to "
                "membrane patch(es) must be tetrahedral and not well-mixed.");

    const auto &comp_tets = comp->_getAllTetIndices();
    tets_set.insert(comp_tets.begin(), comp_tets.end());
  }

  pTet_indices.assign(tets_set.begin(), tets_set.end());
  pTetsN = pTet_indices.size();
  pTrisN = pTri_indices.size();

  ArgErrLogIf(pTrisN == 0, "Membrane contains no triangles.");

  if (verify) {
    verifyMemb();
  }

  // Create sorted set of unique vertex indices
  std::set<vertex_id_t> verts_set;
  for (auto tet : pTet_indices) {
    const auto tet_verts = pTetmesh->_getTet(tet);
    verts_set.insert(tet_verts, tet_verts + 4);
  }

  pVert_indices.assign(verts_set.begin(), verts_set.end());
  pVertsN = pVert_indices.size();

  pTetmesh->_handleMembAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> Memb::isTriInside(const std::vector<index_t> &tri) const {
  return steps::util::map_membership(tri, pTri_indices);
}

////////////////////////////////////////////////////////////////////////////////

void Memb::verifyMemb() {
  CLOG(INFO, "general_log")
      << "Verifying membrane. This can take some time...\n";
  // Now perform some checks on the triangle and print warning if surface is
  // open, or has more than 3 neighbours
  // std::vector<uint>::const_iterator tri_end = pTri_indices.end();
  for (auto const &tri : pTri_indices) {
    uint neighbours = 0;

    // Fetch all triangle neighbours- we don't know how many neighbours there
    // are
    auto const &tritrineighbs = pTetmesh->getTriTriNeighbs(tri);

    for (auto const &t : pTri_indices) {
      // continue of this is the same triangle
      if (t == tri) {
        continue;
      }
      if (tritrineighbs.find(t.get()) != tritrineighbs.end()) {
        neighbours += 1;
      }
    }

    // If 3 neighbours weren't found this is not a closed surface
    // MOST SURFACE TRIANGLES WILL HAVE 3 NEIGHBOURS, BUT
    // IF A TRIANGLE IS ON A SHARP EDGE IT CAN HAVE 5 OR EVEN 7 NEIGHBOURS
    if ((!open()) && neighbours < 3) {
      CLOG(INFO, "general_log") << "WARNING: Triangular surface provided to "
                                   "Membrane initializer function";
      CLOG(INFO, "general_log") << " is thought to be an open surface.\n";
      pOpen = true;
    }
    // What can we say if tri has more than 3 neighbours? Probably a bad
    // surface, though with particularly 'pointed' regions neighbours can be 5,
    // 7, or even higher even number and still be a closed surface. If one of
    // these triangles further happens to be on an edge may have 4, 6, or higher
    // even number, though these numbers would usually suggest a bad surface.
    // However, don't complain yet- perform the multiple surface test later to
    // catch that.
    if (neighbours > 3) {
      CLOG(INFO, "general_log") << "WARNING: Triangle " << tri << " has "
                                << neighbours << " neighbours.\n";
    }
  }

  // Now to perform a basic multiple surfaces test
  {
    // Create the set of triangle neighbours. If this set contains all triangle
    // indices at the end of the loop, this is a single closed surface
    std::set<triangle_id_t> tri_set;
    // Create the queue of triangles- each triangle's neighbours will be
    // added to the queue only once
    std::queue<triangle_id_t> tri_queue;
    // Let's start with the zeroth element
    tri_queue.push(pTri_indices[0]);
    while (!tri_queue.empty()) {
      auto triidx = tri_queue.front();
      auto const &trineighbs = pTetmesh->getTriTriNeighbs(triidx);
      for (auto const &t : pTri_indices) {
        // continue if this is the same triangle
        if (t == triidx) {
          continue;
        }
        if (trineighbs.find(t.get()) != trineighbs.end()) {
          // A triangle's neighbour. Now need to check if it is already in the
          // queue
          if (tri_set.find(t) == tri_set.end()) {
            tri_queue.push(t);
            tri_set.insert(t);
          }
        }
      }
      tri_queue.pop();
    }
    ArgErrLogIf(tri_set.size() != pTri_indices.size(),
                "Triangular surface provided to Membrane initializer function "
                "is a multiple surface.");
  }

  // Loop over all triangles and test that all triangles have one tetrahedron
  // neighbour in the volume.
  for (auto const &tri : pTri_indices) {
    uint tetneighbours = 0;
    // Fetch the two triangle tet neighbours
    const auto *tritetneighbs = pTetmesh->_getTriTetNeighb(tri);

    for (auto t : pTet_indices) {
      if (tritetneighbs[0] == t || tritetneighbs[1] == t) {
        tetneighbours += 1;
      }
    }

    ArgErrLogIf(tetneighbours < 1,
                "Conduction volume provided to Membrane initializer function "
                " is not connected to membrane triangle # " +
                    std::to_string(tri) + ".");
    ArgErrLogIf(tetneighbours > 1,
                "Conduction volume provided to Membrane initializer function "
                " is connected twice to membrane triangle # " +
                    std::to_string(tri) + ".");
  }

  // Now perform the 'multiple volume' test

  // Create the set of tetrahedron's neighbours. If this set contains all
  // tetrahedrons at the end of the loop, this is a single volume.
  std::set<tetrahedron_id_t> tet_set;

  {
    // Create the queue of tetrahedrons- each tetrahedrons' neighbours
    // will be added to the queue only once
    std::queue<tetrahedron_id_t> tet_queue;
    // Start with the zeroth element
    tet_queue.push(pTet_indices[0]);
    while (!tet_queue.empty()) {
      auto tetidx = tet_queue.front();
      const auto *tetneighbs = pTetmesh->_getTetTetNeighb(tetidx);
      for (const auto &t : pTet_indices) {
        // continue if this is the same tet
        if (t == tetidx) {
          continue;
        }
        for (auto i = 0u; i < 4; ++i) {
          if (tetneighbs[i] == t) {
            if (tet_set.insert(t).second) {
              tet_queue.push(t);
            }
            // No need to test other neighbours
            break;
          }
        }
      }
      tet_queue.pop();
    }
  }

  ArgErrLogIf(tet_set.size() != pTet_indices.size(),
              "Conduction volume provided to Membrane initializer function is "
              "a multiple volume.\n");

  // Now to set up the 'virtual membrane triangles'. That is, with an open
  // surface, we need to know which triangles are in the open region.
  if (open()) {
    std::vector<triangle_id_t> trivirttemp;
    // Loop over tetrahedrons and add triangles from surface tets that
    // are not already in pTri_indices
    for (auto const &tet : pTet_indices) {
      // Fetch neighbours
      auto const &tettetneighbs = pTetmesh->getTetTetNeighb(tet);
      // Loop over all tets again and find how many neighbours this tet has.
      // If it has less than 4 it is a candidate for a surface tet whose
      // triangle (in the correct direction) makes up part of the
      // 'virtual surface'
      // std::array<bool, 4> gotneighb{false, false, false, false};
      bool gotneighb[4] = {false, false, false, false};
      for (const auto &tettemp : pTet_indices) {
        if (tettemp == tet) {
          continue;
        }
        for (int i = 0; i < 4; ++i) {
          if (tettetneighbs[i] == tettemp) {
            AssertLog(!gotneighb[i]);
            gotneighb[i] = true;
            break;
          }
        }
      }
      // Now need to check if the triangle in this direction is in
      // pTri_indices
      auto const &tettrineighb = pTetmesh->getTetTriNeighb(tet);
      for (uint j = 0; j < 4; ++j) {
        if (!gotneighb[j]) {
          // Now need to check if the triangle in this direction is in
          // pTri_indices
          bool membtri = false;
          for (auto const &tri : pTri_indices) {
            if (tettrineighb[j] == tri) {
              membtri = true;
              break;
            }
          }
          if (!membtri) {
            trivirttemp.push_back(pTetmesh->getTetTriNeighb(tet)[j]);
          }
        }
      }
    }

    // Ok, let's now go on a little walk and add the correct surface triangles
    // Create the queue of triangles- each triangle's neighbours will be
    // added to the queue only once
    std::queue<triangle_id_t> tri_queue;
    // Let's start with the zeroth element
    tri_queue.push(pTri_indices[0]);

    // Create the set of triangle neighbours. This set will contain all the
    // surface triangles at the end- real and virtual
    std::set<triangle_id_t> tri_set;
    // The virtual triangles, in a set
    std::set<triangle_id_t> trivirt_set;

    while (!tri_queue.empty()) {
      auto triidx = tri_queue.front();
      auto const &trineighbs = pTetmesh->getTriTriNeighbs(triidx);

      for (const auto &t : pTri_indices) {
        // continue if this is the same triangle
        if (t == triidx) {
          continue;
        }
        if (trineighbs.find(t.get()) != trineighbs.end()) {
          // A triangle's neighbour. Now need to check if it is already in the
          // queue
          if (tri_set.find(t) == tri_set.end()) {
            tri_queue.push(t);
            tri_set.insert(t);
          }
        }
      }
      for (const auto &tv : trivirttemp) {
        if (tv == triidx) {
          continue;
        }
        if (trineighbs.find(tv.get()) != trineighbs.end()) {
          if (tri_set.find(tv) == tri_set.end()) {
            tri_queue.push(tv);
            tri_set.insert(tv);
            trivirt_set.insert(tv);
          }
        }
      }
      tri_queue.pop();
    }

    // So trivirt_set holds all the virtual triangles. Just need to copy to
    // vector
    for (const auto &tvset : trivirt_set) {
      pTrivirt_indices.push_back(tvset);
    }

    pTriVirtsN = pTrivirt_indices.size();
  }
}

} // namespace tetmesh
} // namespace steps

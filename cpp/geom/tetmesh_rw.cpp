////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */


// STL headers.
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "comp.hpp"
#include "patch.hpp"
#include "tetmesh_rw.hpp"
#include "tetmesh.hpp"
#include "tmcomp.hpp"
#include "tmpatch.hpp"

USING(std, endl);
USING(std, map);
USING(std, ifstream);
USING(std, ofstream);
USING(std, ostringstream);
USING(std, set);
USING(std, string);
USING(std, vector);
USING(steps::tetmesh, Tetmesh);
USING(steps::tetmesh, TmComp);
USING(steps::tetmesh, TmPatch);
USING(steps::wm, Comp);
USING(steps::wm, Patch);

////////////////////////////////////////////////////////////////////////////////

// TODO:
// * Change to a binary format, this will obviously be much more effective
//   for larger meshes. Maybe use something like the NetCDF library for
//   this. It will make loading the mesh from e.g. Matlab easier.
// * Add a LOT more error checking and throw exceptions (while keeping in
//   mind to clean up along the call path!! As always!!!!!!!!)

////////////////////////////////////////////////////////////////////////////////

Tetmesh * steps::tetmesh::loadASCII(string pathname)
{
    typedef vector<string> strvec_t;
    typedef strvec_t::const_iterator strvec_ci;
    typedef vector<double> doublevec_t;
    typedef vector<uint> uintvec_t;
    typedef map<string, TmComp*> compmap_t;

    ifstream mf(pathname.c_str());
    if (!mf)
    {
        ostringstream os;
        os << "Cannot open file \"" << pathname << "\"";
        throw steps::IOErr(os.str());
    }

    // Read vertices.
    uint nverts = 0;
    mf >> nverts;
    doublevec_t verts(nverts * 3);
    for (uint i = 0; i < nverts * 3; ++i)
    {
        double v;
        mf >> v;
        verts[i] = v;
    }

    // Read triangles.
    uint ntris = 0;
    mf >> ntris;
    uintvec_t tris(ntris * 3);
    for (uint i = 0; i < ntris * 3; ++i)
    {
        uint v;
        mf >> v;
        tris[i] = v;
    }

    // Read tetrahedrons.
    uint ntets = 0;
    mf >> ntets;
    uintvec_t tets(ntets * 4);
    for (uint i = 0; i < ntets * 4; ++i)
    {
        uint v;
        mf >> v;
        tets[i] = v;
    }

    // Build mesh and set stuff using new style constructor).
    Tetmesh * m = new Tetmesh(verts, tets, tris);

    // Read compartments.
    uint ncomps = 0;
    mf >> ncomps;
    compmap_t compmap;
    for (uint c = 0; c < ncomps; ++c)
    {
        string compid;
        mf >> compid;

        uint nvolsys;
        mf >> nvolsys;
        strvec_t volsys;
        for (uint v = 0; v < nvolsys; ++v)
        {
            string volsysid;
            mf >> volsysid;
            volsys.push_back(volsysid);
        }

        uint ntets_in_c;
        mf >> ntets_in_c;
        uintvec_t comptets(ntets_in_c);
        for (uint t = 0; t < ntets_in_c; ++t)
        {
            uint v;
            mf >> v;
            comptets[t] = v;
        }

        TmComp * newcomp = new TmComp(compid, m, comptets);
        compmap[compid] = newcomp;
        strvec_ci volsys_end = volsys.end();
        for (strvec_ci vsys = volsys.begin(); vsys != volsys_end; ++vsys)
            newcomp->addVolsys(*vsys);
    }

    // Read patches.
    uint npatches = 0;
    mf >> npatches;
    for (uint p = 0; p < npatches; ++p)
    {
        // Read patch name.
        string patchid;
        mf >> patchid;

        // Read inner comp.
        uint icomp_switch = 0;
        TmComp * icomp = 0;
        mf >> icomp_switch;
        if (icomp_switch != 0)
        {
            string icomp_id;
            mf >> icomp_id;
            icomp = compmap[icomp_id];
        }

        // Read outer comp.
        uint ocomp_switch = 0;
        TmComp * ocomp = 0;
        mf >> ocomp_switch;
        if (ocomp_switch != 0)
        {
            string ocomp_id;
            mf >> ocomp_id;
            ocomp = compmap[ocomp_id];
        }

        // Read surface systems.
        uint nsurfsys = 0;
        mf >> nsurfsys;
        strvec_t surfsys;
        for (uint s = 0; s < nsurfsys; ++s)
        {
            string ssysid;
            mf >> ssysid;
            surfsys.push_back(ssysid);
        }

        // Read triangles.
        uint ntris_in_p = 0;
        mf >> ntris_in_p;
        uintvec_t patchtris(ntris_in_p);
        for (uint t = 0; t < ntris_in_p; ++t)
        {
            uint v;
            mf >> v;
            patchtris[t] = v;
        }

        TmPatch * patch = new TmPatch(patchid, m, patchtris, icomp, ocomp);
        strvec_ci surfsys_end = surfsys.end();
        for (strvec_ci ssys = surfsys.begin(); ssys != surfsys_end; ++ssys)
            patch->addSurfsys(*ssys);
    }

    mf.close();
    return m;
}

////////////////////////////////////////////////////////////////////////////////

void steps::tetmesh::saveASCII(string pathname, Tetmesh * m)
{
    typedef set<string> strset;
    typedef strset::const_iterator strset_ci;
    typedef vector<uint> uintvec;
    typedef uintvec::const_iterator uintvec_ci;

    if (m == 0)
    {
        ostringstream os;
        os << "No model specified";
        throw steps::ArgErr(os.str());
    }

    ofstream mf(pathname.c_str());
    if (!mf)
    {
        ostringstream os;
        os << "Cannot open file \"" << pathname << "\"";
        throw steps::IOErr(os.str());
    }

    // Increase digit precision a bit to accurately store doubles
    // in ASCII. Otherwise last few binary digits might get rounded
    // wrongly when reading back in.
    mf.precision(13);

    // Dump vertices.
    uint nverts = m->countVertices();
    mf << nverts << endl;
    for (uint i = 0; i < nverts; ++i)
    {
        double * verts = m->_getVertex(i);
        mf.width(20);
        mf << verts[0] << "    ";
        mf.width(20);
        mf << verts[1] << "    ";
        mf.width(20);
        mf << verts[2] << endl;
    }
    mf << endl;

    // Dump triangles.
    uint ntris = m->countTris();
    mf << ntris << endl;
    for (uint i = 0; i < ntris; ++i)
    {
        uint * tri = m->_getTri(i);
        mf.width(8);
        mf << tri[0] << "  ";
        mf.width(8);
        mf << tri[1] << "  ";
        mf.width(8);
        mf << tri[2] << endl;
    }
    mf << endl;

    // Dump tetrahedrons.
    uint ntets = m->countTets();
    mf << ntets << endl;
    for (uint i = 0; i < ntets; ++i)
    {
        uint * tet = m->_getTet(i);
        mf.width(8);
        mf << tet[0] << "  ";
        mf.width(8);
        mf << tet[1] << "  ";
        mf.width(8);
        mf << tet[2] << "  ";
        mf.width(8);
        mf << tet[3] << endl;
    }
    mf << endl;

    // Dump compartments.
    uint ncomps = m->_countComps();
    mf << ncomps << endl;
    for (uint cidx = 0; cidx < ncomps; ++cidx)
    {
        TmComp * comp = m->getTetComp(cidx);
        mf << comp->getID() << endl;

        strset volsys = comp->getVolsys();
        mf << volsys.size() << endl;
        strset_ci v_end = volsys.end();
        for (strset_ci v = volsys.begin(); v != v_end; ++v)
        {
            mf << *v << endl;
        }

        uintvec tets = comp->_getAllTetIndices();
        mf << tets.size() << endl;
        uintvec_ci t_end = tets.end();
        uint numctr = 0;
        for (uintvec_ci t = tets.begin(); t != t_end; ++t)
        {
            mf.width(8);
            mf << *t;
            numctr++;
            if (numctr == 8)
            {
                numctr = 0;
                mf << endl;
            }
            else
            {
                mf << "  ";
            }
        }
        if (numctr != 0) mf << endl;
        mf << endl;
    }

    // Dump patches.
    uint npatches = m->_countPatches();
    mf << npatches << endl;
    for (uint pidx = 0; pidx < npatches; ++pidx)
    {
        TmPatch * patch = m->getTriPatch(pidx);
        mf << patch->getID() << endl;

        Comp * icomp = patch->getIComp();
        if (icomp == 0) mf << "0" << endl;
        else mf << "1" << "    " << icomp->getID() << endl;

        Comp * ocomp = patch->getOComp();
        if (ocomp == 0) mf << "0" << endl;
        else mf << "1" << "    " << ocomp->getID() << endl;

        strset surfsys = patch->getSurfsys();
        mf << surfsys.size() << endl;
        strset_ci s_end = surfsys.end();
        for (strset_ci s = surfsys.begin(); s != s_end; ++s)
        {
            mf << *s << endl;
        }

        uintvec tris = patch->_getAllTriIndices();
        mf << tris.size() << endl;
        uintvec_ci t_end = tris.end();
        uint numctr = 0;
        for (uintvec_ci t = tris.begin(); t != t_end; ++t)
        {
            mf.width(8);
            mf << *t;
            numctr++;
            if (numctr == 8)
            {
                numctr = 0;
                mf << endl;
            }
            else
            {
                mf << "  ";
            }
        }
        if (numctr != 0) mf << endl;
        mf << endl;
    }

    mf.close();
}

////////////////////////////////////////////////////////////////////////////////

// END

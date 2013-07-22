/*!
 * \file solution_template.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author J.B. Scoggins (von Karman Institute for Fluid Dynmaics)
 * \version 1.0.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/solution_structure.hpp"

#ifdef CHECK

CRTESolution::CRTESolution(void) : CSolution() { }

CRTESolution::CRTESolution(CGeometry *geometry, CConfig *config) : CSolution() { }

CRTESolution::~CRTESolution(void) { }

void CRTESolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CNumerics **solver, CConfig *config, unsigned short iRKStep) { }

void CRTESolution::SetTime_Step(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CRTESolution::Centered_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
                                          CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

void CRTESolution::Upwind_Residual(
    CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
    CConfig *config, unsigned short iMesh)
{
    // Loop over each face
    for(int iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		// Points in edge and normal vectors
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		solver->SetNormal(geometry->edge[iEdge]->GetNormal());
        
		// Set conservative variables w/o reconstruction
		solver->SetConservative(
            node[iPoint]->GetSolution(),
            node[jPoint]->GetSolution());
        
        // Compute the residual and add to left cell, subtract from right
        solver->SetResidual(Res_Conv, NULL, NULL, config);
        node[iPoint]->AddRes_Conv(Res_Conv);
		node[jPoint]->SubtractRes_Conv(Res_Conv);
    }
}

void CRTESolution::Source_Residual(
    CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
    CConfig *config, unsigned short iMesh)
{
    // Loop over every point in the domain
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        // Set solution and control volume
        solver->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
        solver->SetVolume(geometry->node[iPoint]->GetVolume());
        
        // Compute Residual
        solver->SetResidual(Residual, NULL, config);
        
        /*--- Subtract Residual ---*/
        node[iPoint]->SubtractRes_Conv(Residual);
    }
}

void CRTESolution::Solve_LinearSystem(CGeometry *geometry, CSolution **solution_container, CConfig *config,
                                           unsigned short iMesh) { }

void CRTESolution::BC_Euler_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
                                      unsigned short val_marker) { }

void CRTESolution::BC_NS_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) { }

void CRTESolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
                                     unsigned short val_marker) { }

void CRTESolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
                                 unsigned short val_marker) { }

void CRTESolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
                                  unsigned short val_marker) { }

void CRTESolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config,
                                     unsigned short val_marker) { }

void CRTESolution::BC_Custom(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short val_marker) { }

void CRTESolution::ExplicitRK_Iteration(CGeometry *geometry, CSolution **solution_container,
                                             CConfig *config, unsigned short iRKStep) { }

void CRTESolution::ExplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

void CRTESolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) { }

#endif

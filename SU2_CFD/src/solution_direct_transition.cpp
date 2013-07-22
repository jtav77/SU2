/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
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

CTransLMSolution::CTransLMSolution(void) : CTurbSolution() {}

CTransLMSolution::CTransLMSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolution() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double Density_Inf, Viscosity_Inf, tu_Inf, nu_tilde_Inf, Factor_nu_Inf, dull_val, rey, mach;
	ifstream restart_file;
	char *cstr;
	string text_line;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem --> 2 Transport equations (intermittency,Reth) ---*/
	nVar = 2;
	
	if (iMesh == MESH_0) {
		
		/*--- Define some auxillary vectors related to the residual ---*/
		Residual     = new double[nVar]; Residual_RMS = new double[nVar];
		Residual_i   = new double[nVar]; Residual_j   = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max    = new unsigned long[nVar];
    
    
		/*--- Define some auxiliar vector related with the solution ---*/
		Solution   = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];
		
		/*--- Define some auxiliar vector related with the flow solution ---*/
		FlowSolution_i = new double [nDim+2]; FlowSolution_j = new double [nDim+2];
		
    xsol = new double [geometry->GetnPoint()*nVar];
    xres = new double [geometry->GetnPoint()*nVar];
    
		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {
      
			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar];
			Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar];
				Jacobian_j[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Initialize_SparseMatrix_Structure(&Jacobian, nVar, nVar, geometry, config);
      
		}
    
    /*--- Computation of gradients by least squares ---*/
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new double [nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
      cvector = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        cvector[iVar] = new double [nDim];
    }
    
    /*--- Read farfield conditions from config ---*/
    Density_Inf       = config->GetDensity_FreeStreamND();
    Viscosity_Inf     = config->GetViscosity_FreeStreamND();
    Intermittency_Inf = config->GetIntermittency_FreeStream();
    tu_Inf            = config->GetTurbulenceIntensity_FreeStream();
    
    /*-- Initialize REth from correlation --*/
    if (tu_Inf <= 1.3) {
      REth_Inf = (1173.51-589.428*tu_Inf+0.2196/(tu_Inf*tu_Inf));
    } else {
      REth_Inf = 331.5*pow(tu_Inf-0.5658,-0.671);
    }
    rey = config->GetReynolds();
    mach = config->GetMach_FreeStreamND();
    
    //  REth_Inf *= mach/rey;
    cout << "REth_Inf = " << REth_Inf << ", rey: "<< rey << " -AA" << endl;
    
    /*--- Factor_nu_Inf in [3.0, 5.0] ---*/
    Factor_nu_Inf = config->GetNuFactor_FreeStream();
    nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
    
    /*--- Eddy viscosity ---*/
    double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
    double muT_Inf;
    Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
    Ji_3 = Ji*Ji*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    muT_Inf = Density_Inf*fv1*nu_tilde_Inf;
    
    /*--- Restart the solution from file information ---*/
    if (!restart) {
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        // TODO: Erase this bubble of specially initialized points -AA
        node[iPoint] = new CTransLMVariable(nu_tilde_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);
      }
    }
    
	}
	else {
    cout << "No LM restart yet!!" << endl; // TODO, Aniket
    int j;
    cin >> j;
		string mesh_filename = config->GetSolution_FlowFileName();
		cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			getline(restart_file,text_line);
			istringstream point_line(text_line);
			if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
			node[iPoint] = new CTurbSAVariable(Solution[0], 0, nDim, nVar, config);
		}
		restart_file.close();
	}
  
}

CTransLMSolution::~CTransLMSolution(void){
	unsigned short iVar, iDim;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Solution;
	delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i; delete [] Vector_j;
	delete [] FlowSolution_i; delete [] FlowSolution_j;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i; delete [] Jacobian_j;
	
	delete [] xsol; delete [] xres;
	
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
}

void CTransLMSolution::Preprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		Set_Residual_Zero(iPoint);
  Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CTransLMSolution::Postprocessing(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iMesh) {

  /*--- Correction for separation-induced transition, Replace intermittency with gamma_eff ---*/
  // Moved this to Source_Residual()
	//for (unsigned int iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
  //  node[iPoint]->SetGammaEff();
}

void CTransLMSolution::ImplicitEuler_Iteration(CGeometry *geometry, CSolution **solution_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, Vol;
  
  
  /*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    Delta_flow = Vol / (solution_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Delta = Delta_flow;
    Jacobian.AddVal2Diag(iPoint,Delta);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      
      /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
      xres[total_index] = -xres[total_index];
      xsol[total_index] = 0.0;
      AddRes_RMS(iVar, xres[total_index]*xres[total_index]*Vol);
      AddRes_Max(iVar, fabs(xres[total_index]), geometry->node[iPoint]->GetGlobalIndex());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      xres[total_index] = 0.0;
      xsol[total_index] = 0.0;
    }
  }
  
	/*--- Solve the linear system (Stationary iterative methods) ---*/
	if (config->GetKind_Linear_Solver() == SYM_GAUSS_SEIDEL)
		Jacobian.SGSSolution(xres, xsol, config->GetLinear_Solver_Error(),
                         config->GetLinear_Solver_Iter(), false, geometry, config);
	
	if (config->GetKind_Linear_Solver() == LU_SGS)
    Jacobian.LU_SGSIteration(xres, xsol, geometry, config);
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
	if ((config->GetKind_Linear_Solver() == BCGSTAB) ||
      (config->GetKind_Linear_Solver() == GMRES)) {
		
		CSysVector rhs_vec((const unsigned int)geometry->GetnPoint(),
                       (const unsigned int)geometry->GetnPointDomain(), nVar, xres);
		CSysVector sol_vec((const unsigned int)geometry->GetnPoint(),
                       (const unsigned int)geometry->GetnPointDomain(), nVar, xsol);
		
		CMatrixVectorProduct* mat_vec = new CSparseMatrixVectorProduct(Jacobian, geometry, config);
    
		CPreconditioner* precond = NULL;
		if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
			Jacobian.BuildJacobiPreconditioner();
			precond = new CJacobiPreconditioner(Jacobian, geometry, config);
		}
		else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
			Jacobian.BuildJacobiPreconditioner();
			precond = new CLineletPreconditioner(Jacobian, geometry, config);
		}
		else if (config->GetKind_Linear_Solver_Prec() == NO_PREC) {
			precond = new CIdentityPreconditioner(Jacobian, geometry, config);
		}
    
		CSysSolve system;
		if (config->GetKind_Linear_Solver() == BCGSTAB)
			system.BCGSTAB(rhs_vec, sol_vec, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                     config->GetLinear_Solver_Iter(), false);
		else if (config->GetKind_Linear_Solver() == GMRES)
			system.GMRES(rhs_vec, sol_vec, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
		
		sol_vec.CopyToArray(xsol);
		delete mat_vec;
		delete precond;
    
	}
	
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++)
      node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*xsol[iPoint*nVar+iVar]);
  }
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CTransLMSolution::Upwind_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CConfig *config, unsigned short iMesh) {
  double *trans_var_i, *trans_var_j, *U_i, *U_j;
  unsigned long iEdge, iPoint, jPoint;

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge and normal vectors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    solver->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solution_container[FLOW_SOL]->node[iPoint]->GetSolution();
    U_j = solution_container[FLOW_SOL]->node[jPoint]->GetSolution();
    solver->SetConservative(U_i, U_j);

    /*--- Transition variables w/o reconstruction ---*/
    trans_var_i = node[iPoint]->GetSolution();
    trans_var_j = node[jPoint]->GetSolution();
    solver->SetTransVar(trans_var_i, trans_var_j);

    /*--- Add and subtract Residual ---*/
    solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    AddResidual(iPoint, Residual);
    SubtractResidual(jPoint, Residual);

    /*--- Implicit part ---*/
    Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
    Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
    Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
    Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);

  }

}


void CTransLMSolution::Viscous_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																				CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    solver->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());
    
    solver->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                            solution_container[FLOW_SOL]->node[jPoint]->GetSolution());
    
    /*--- Laminar Viscosity ---*/
    solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                solution_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    /*--- Eddy Viscosity ---*/
    solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                             solution_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());
    
    /*--- Transition variables w/o reconstruction, and its gradients ---*/
    solver->SetTransVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    solver->SetTransVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    solver->SetConsVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient(),
                               solution_container[FLOW_SOL]->node[jPoint]->GetGradient());
    
    
    /*--- Compute residual, and Jacobians ---*/
    solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    SubtractResidual(iPoint, Residual);
    AddResidual(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTransLMSolution::Source_Residual(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *second_solver,
																			 CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  double gamma_sep;

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		
		/*--- Conservative variables w/o reconstruction ---*/
		solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
		
		/*--- Gradient of the primitive and conservative variables ---*/
		solver->SetPrimVarGradient(solution_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		
		/*--- Laminar and eddy viscosity ---*/
		solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
    solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);
		
		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		solver->SetTransVar(node[iPoint]->GetSolution(), NULL);
		
		/*--- Set volume ---*/
		solver->SetVolume(geometry->node[iPoint]->GetVolume());
		
		/*--- Set distance to the surface ---*/
		solver->SetDistance(geometry->node[iPoint]->GetWallDistance(), 0.0);
		
		/*--- Compute the source term ---*/
		solver->SetResidual_TransLM(Residual, Jacobian_i, NULL, config, gamma_sep);
		
    /*-- Store gamma_sep in variable class --*/
    node[iPoint]->SetGammaSep(gamma_sep);
    node[iPoint]->SetGammaEff();

		/*--- Subtract residual and the jacobian ---*/
		SubtractResidual(iPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

	}
  
}

void CTransLMSolution::Source_Template(CGeometry *geometry, CSolution **solution_container, CNumerics *solver,
																			 CConfig *config, unsigned short iMesh) {
}

void CTransLMSolution::BC_HeatFlux_Wall(CGeometry *geometry, CSolution **solution_container, CNumerics *solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	unsigned short iVar;
  int total_index;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

      /* --- Impose boundary values (Dirichlet) ---*/
      Solution[0] = 0.0;
      Solution[1] = 0.0;
			node[iPoint]->SetSolution_Old(Solution);
			Set_Residual_Zero(iPoint);

			/*--- includes 1 in the diagonal ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
		}
	}
  
}

void CTransLMSolution::BC_Far_Field(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) { }

void CTransLMSolution::BC_Inlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();
	unsigned long iVertex, iPoint, Point_Normal;
	double P_Total, T_Total, Velocity[3];
	double Velocity2, H_Total, Temperature, Riemann;
	double Pressure, Density, Energy, *Flow_Dir, Mach2;
	double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
	double alpha, aa, bb, cc, dd;
	double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
	double Gas_Constant = config->GetGas_ConstantND();
	double *Normal = new double[nDim];
  
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
	bool incompressible = config->GetIncompressible();
	string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
			double Area = 0.0; double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;
      
			/*--- Current conservative variables at this boundary node (U_domain) ---*/
			for (iVar = 0; iVar < solution_container[FLOW_SOL]->GetnVar(); iVar++)
				FlowSolution_i[iVar] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);
      
			/*--- Build the fictitious intlet state based on characteristics ---*/
			if (incompressible) {
        
				/*--- Pressure computation using the internal value ---*/
				FlowSolution_j[0] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
        
				/*--- The velocity is computed from the interior, and normal
         derivative for the density ---*/
				for (iDim = 0; iDim < nDim; iDim++)
					FlowSolution_j[iDim+1] = solution_container[FLOW_SOL]->GetVelocity_Inf(iDim)*solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc();
        
			}
			else {
        
				/*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/
        
				switch (Kind_Inlet) {
            
            /*--- Total properties have been specified at the inlet. ---*/
          case TOTAL_CONDITIONS:
            
            /*--- Retrieve the specified total conditions for this inlet. ---*/
            P_Total  = config->GetInlet_Ptotal(Marker_Tag);
            T_Total  = config->GetInlet_Ttotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
            
            /*--- Non-dim. the inputs if necessary. ---*/
            P_Total /= config->GetPressure_Ref();
            T_Total /= config->GetTemperature_Ref();
            
            /*--- Store primitives and set some variables for clarity. ---*/
            Density = FlowSolution_i[0];
            Velocity2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
              Velocity2 += Velocity[iDim]*Velocity[iDim];
            }
            Energy      = FlowSolution_i[nVar-1]/Density;
            Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
            H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
            SoundSpeed2 = Gamma*Pressure/Density;
            
            /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
            Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
            for (iDim = 0; iDim < nDim; iDim++)
              Riemann += Velocity[iDim]*UnitaryNormal[iDim];
            
            /*--- Total speed of sound ---*/
            SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy
                                                            + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
            
            /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/
            alpha = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];
            
            /*--- Coefficients in the quadratic equation for the velocity ---*/
            aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
            bb = -1.0*Gamma_Minus_One*alpha*Riemann;
            cc =  0.5*Gamma_Minus_One*Riemann*Riemann
            -2.0*SoundSpeed_Total2/Gamma_Minus_One;
            
            /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/
            dd = bb*bb - 4.0*aa*cc;
            dd = sqrt(max(0.0,dd));
            Vel_Mag   = (-bb + dd)/(2.0*aa);
            Vel_Mag   = max(0.0,Vel_Mag);
            Velocity2 = Vel_Mag*Vel_Mag;
            
            /*--- Compute speed of sound from total speed of sound eqn. ---*/
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
            
            /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
            Mach2 = Velocity2/SoundSpeed2;
            Mach2 = min(1.0,Mach2);
            Velocity2   = Mach2*SoundSpeed2;
            Vel_Mag     = sqrt(Velocity2);
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
            
            /*--- Compute new velocity vector at the inlet ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
            
            /*--- Static temperature from the speed of sound relation ---*/
            Temperature = SoundSpeed2/(Gamma*Gas_Constant);
            
            /*--- Static pressure using isentropic relation at a point ---*/
            Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);
            
            /*--- Density at the inlet from the gas law ---*/
            Density = Pressure/(Gas_Constant*Temperature);
            
            /*--- Using pressure, density, & velocity, compute the energy ---*/
            Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
            
            /*--- Conservative variables, using the derived quantities ---*/
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Velocity[0]*Density;
            FlowSolution_j[2] = Velocity[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Velocity[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }
            
            break;
            
            /*--- Mass flow has been specified at the inlet. ---*/
          case MASS_FLOW:
            
            /*--- Retrieve the specified mass flow for the inlet. ---*/
            Density  = config->GetInlet_Ttotal(Marker_Tag);
            Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
            
            /*--- Non-dim. the inputs if necessary. ---*/
            Density /= config->GetDensity_Ref();
            Vel_Mag /= config->GetVelocity_Ref();
            
            /*--- Get primitives from current inlet state. ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = node[iPoint]->GetVelocity(iDim, incompressible);
            Pressure    = node[iPoint]->GetPressure(incompressible);
            SoundSpeed2 = Gamma*Pressure/FlowSolution_i[0];
            
            /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
            Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
            for (iDim = 0; iDim < nDim; iDim++)
              Riemann += Velocity[iDim]*UnitaryNormal[iDim];
            
            /*--- Speed of sound squared for fictitious inlet state ---*/
            SoundSpeed2 = Riemann;
            for (iDim = 0; iDim < nDim; iDim++)
              SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitaryNormal[iDim];
            
            SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
            SoundSpeed2 = SoundSpeed2*SoundSpeed2;
            
            /*--- Pressure for the fictitious inlet state ---*/
            Pressure = SoundSpeed2*Density/Gamma;
            
            /*--- Energy for the fictitious inlet state ---*/
            Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Vel_Mag*Vel_Mag;
            
            /*--- Conservative variables, using the derived quantities ---*/
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Vel_Mag*Flow_Dir[0]*Density;
            FlowSolution_j[2] = Vel_Mag*Flow_Dir[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Vel_Mag*Flow_Dir[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }
            
            break;
				}
			}
      
			/*--- Set the conservative variable states. Note that we really only
       need the density and momentum components so that we can compute
       the velocity for the convective part of the upwind residual. ---*/
			conv_solver->SetConservative(FlowSolution_i, FlowSolution_j);
      
			/*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = Intermittency_Inf;

      Solution_i[1] = node[iPoint]->GetSolution(1);
      Solution_j[1] = REth_Inf;
      
			conv_solver->SetTransVar(Solution_i, Solution_j);
      
			/*--- Set various other quantities in the conv_solver class ---*/
			conv_solver->SetNormal(Normal);
      
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(),
                               geometry->node[iPoint]->GetRotVel());
				conv_solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
			/*--- Compute the residual using an upwind scheme ---*/
      conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
      AddResidual(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution ---*/
//			visc_solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//			visc_solver->SetNormal(Normal);
//      
//			/*--- Conservative variables w/o reconstruction ---*/
//			visc_solver->SetConservative(FlowSolution_i, FlowSolution_j);
//      
//			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
//			if (incompressible) {
//				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
//                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
//				visc_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
//                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
//			} else {
//				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
//                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
//			}
//      
//			/*--- Eddy Viscosity ---*/
//			visc_solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
//                                    solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
//      
//			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//			visc_solver->SetTurbVar(Solution_i, Solution_j);
//			visc_solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//      
//			/*--- Compute residual, and Jacobians ---*/
//			visc_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//      
//			/*--- Subtract residual, and update Jacobians ---*/
//			SubtractResidual(iPoint, Residual);
//			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
  
}

void CTransLMSolution::BC_Outlet(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver,
																 CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
  
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool grid_movement  = config->GetGrid_Movement();
  
	double *Normal = new double[nDim];
  
	/*--- Loop over all the vertices on this boundary marker ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
		/*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
			Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
			/*--- Set the conservative variables, density & Velocity same as in
			 the interior, don't need to modify pressure for the convection. ---*/
			conv_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
			/*--- Set the turbulent variables. Here we use a Neumann BC such
			 that the turbulent variable is copied from the interior of the
			 domain to the outlet before computing the residual.
			 Solution_i --> TurbVar_internal,
			 Solution_j --> TurbVar_outlet ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
				Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
			}
			conv_solver->SetTransVar(Solution_i, Solution_j);
      
			/*--- Set Normal (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++)
				Normal[iDim] = -Normal[iDim];
			conv_solver->SetNormal(Normal);
      
			/*--- Set various quantities in the solver class ---*/
			if (incompressible)
				conv_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
			if (rotating_frame) {
				conv_solver->SetRotVel(geometry->node[iPoint]->GetRotVel(),
                               geometry->node[iPoint]->GetRotVel());
				conv_solver->SetRotFlux(-geometry->vertex[val_marker][iVertex]->GetRotFlux());
			}
			if (grid_movement)
				conv_solver->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
			/*--- Compute the residual using an upwind scheme ---*/
			conv_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
			AddResidual(iPoint, Residual);
      
			/*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution ---*/
//			visc_solver->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//			visc_solver->SetNormal(Normal);
//      
//			/*--- Conservative variables w/o reconstruction ---*/
//			visc_solver->SetConservative(solution_container[FLOW_SOL]->node[iPoint]->GetSolution(),
//                                   solution_container[FLOW_SOL]->node[iPoint]->GetSolution());
//      
//			/*--- Laminar Viscosity, and density (incompresible solver)  ---*/
//			if (incompressible) {
//				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
//                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
//				visc_solver->SetDensityInc(solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
//                                   solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
//			} else {
//				visc_solver->SetLaminarViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
//                                         solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
//			}
//      
//			/*--- Eddy Viscosity ---*/
//			visc_solver->SetEddyViscosity(solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
//                                    solution_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
//      
//			/*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//			visc_solver->SetTurbVar(Solution_i, Solution_j);
//			visc_solver->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//      
//			/*--- Compute residual, and Jacobians ---*/
//			visc_solver->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//      
//			/*--- Subtract residual, and update Jacobians ---*/
//			SubtractResidual(iPoint, Residual);
//			Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
		}
	}
  
	/*--- Free locally allocated memory ---*/
	delete[] Normal;
  
}

void CTransLMSolution::BC_Sym_Plane(CGeometry *geometry, CSolution **solution_container, CNumerics *conv_solver, CNumerics *visc_solver, CConfig *config, unsigned short val_marker) {
	/*--- Convective fluxes across symmetry plane are equal to zero. ---*/
}

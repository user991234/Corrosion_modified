// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
variableAttributeLoader::loadVariableAttributes()
{
  // ---------------------------
  // Original 8 corrosion vars
  // ---------------------------

  // Variable 0
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n, grad(psi), irxn");
  set_dependencies_gradient_term_RHS(0, "psi, grad(mu), irxn");

  // Variable 1
  set_variable_name(1, "mu");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "n, psi");
  set_dependencies_gradient_term_RHS(1, "grad(n)");

  // Variable 2
  set_variable_name(2, "psi");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "psi, grad(psi), irxn");
  set_dependencies_gradient_term_RHS(2, "psi, grad(mupsi), irxn");

  // Variable 3
  set_variable_name(3, "mupsi");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, AUXILIARY);

  set_dependencies_value_term_RHS(3, "n, psi");
  set_dependencies_gradient_term_RHS(3, "grad(psi)");

  // Variable 4
  set_variable_name(4, "cM");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(4, "cM, grad(cM), psi, grad(psi), irxn");
  set_dependencies_gradient_term_RHS(4, "cM, grad(cM), grad(Phi)");

  // Variable 5
  set_variable_name(5, "cP");
  set_variable_type(5, SCALAR);
  set_variable_equation_type(5, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(5, "cP, grad(cP), psi, grad(psi)");
  set_dependencies_gradient_term_RHS(5, "cP, grad(cP), grad(Phi)");

  // Variable 6
  set_variable_name(6, "Phi");
  set_variable_type(6, SCALAR);
  set_variable_equation_type(6, TIME_INDEPENDENT);

  set_dependencies_value_term_LHS(6,
    "n, psi, grad(psi), cM, grad(cM), cP, Phi, grad(Phi), change(Phi), irxn");
  set_dependencies_gradient_term_LHS(6, "n, psi, cM, cP, grad(change(Phi))");
  set_dependencies_value_term_RHS(6, "grad(psi), irxn");
  set_dependencies_gradient_term_RHS(6, "psi, grad(Phi), grad(cM), grad(cP)");

  // Variable 7
  set_variable_name(7, "irxn");
  set_variable_type(7, SCALAR);
  set_variable_equation_type(7, AUXILIARY);

  set_dependencies_value_term_RHS(7, "cM, cP, Phi");
  set_dependencies_gradient_term_RHS(7, "");

  // ------------------------------------------
  // NEW: 8 additional grain fields g0..g7
  // ------------------------------------------
  // We'll label them indices 8..15
  // Each is "EXPLICIT_TIME_DEPENDENT" with typical Allen–Cahn references:
  for (unsigned int i = 8; i < 16; i++)
    {
      unsigned int grain_index = i - 8;  // 0..7
      std::string var_name = "g" + std::to_string(grain_index);

      set_variable_name(i, var_name);
      set_variable_type(i, SCALAR);
      set_variable_equation_type(i, EXPLICIT_TIME_DEPENDENT);

      // We might need the values of all g0..g7 if they share cross-coupling:
      set_dependencies_value_term_RHS(i, "g0, g1, g2, g3, g4, g5, g6, g7");
      // We might also need the gradients for K-term (Laplacian):
      set_dependencies_gradient_term_RHS(i,
        "grad(g0), grad(g1), grad(g2), grad(g3), grad(g4), grad(g5), grad(g6), grad(g7)");
    }
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time-dependent)
// =============================================================================================
// This function calculates the right-hand side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index
// for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // ---------------------------
  // (A) Original corrosion code
  // ---------------------------
  // The following lines are the existing logic for variables n, psi, cM, cP, etc.
  // We do NOT modify them; we just keep them as-is.

  // Timestep
  scalarvalueType delt = constV(userInputs.dtValue);

  // The order parameter and its derivatives
  scalarvalueType n = variable_list.get_scalar_value(0);

  // The chemical potential and its derivatives
  scalargradType mux = variable_list.get_scalar_gradient(1);

  // The domain parameter and its derivatives
  scalarvalueType psi  = variable_list.get_scalar_value(2);
  scalargradType  psix = variable_list.get_scalar_gradient(2);

  // The chemical potential of the domain parameter
  scalargradType mupsix = variable_list.get_scalar_gradient(3);

  // The concentration of metal ion and its derivatives
  scalarvalueType cM  = variable_list.get_scalar_value(4);
  scalargradType  cMx = variable_list.get_scalar_gradient(4);

  // The concentration of supporting cation and its derivatives
  scalarvalueType cP  = variable_list.get_scalar_value(5);
  scalargradType  cPx = variable_list.get_scalar_gradient(5);

  // The electrolite potential
  scalargradType Phix = variable_list.get_scalar_gradient(6);

  // The reaction current
  scalarvalueType irxn = variable_list.get_scalar_value(7);

  // (Same capping logic as before)
  psi = std::min(psi, constV(1.0));
  psi = std::max(psi, constV(lthresh));
  n   = std::min(n, constV(1.0));
  n   = std::max(n, constV(lthresh));

  // (Same calculations for the corrosion PDEs as shown in your original code)
  // ...
  // Submitting the terms for n, psi, cM, cP:
  //    variable_list.set_scalar_value_term_RHS(0, rnV);
  //    variable_list.set_scalar_gradient_term_RHS(0, rnxV);
  //    variable_list.set_scalar_value_term_RHS(2, rpsiV);
  //    ...
  //    etc.

  // ----------------------------
  // (B) New logic for g0..g7
  // ----------------------------
  // Indices 8..15 => g0..g7
  // We'll do a minimal Allen–Cahn-like update:
  //   dg_i/dt = - Mgrain [ -W*g_i + W*g_i^3 + ... - Kgrain Laplacian(g_i)]
  // For demonstration, let's do a simple double-well for each g_i with no cross coupling:
  //   dF/dg_i = -W*g_i + W*g_i^3 - Kgrain \nabla^2 g_i
  // You can add cross terms if needed.

  // First, gather the scalar values and gradients for g0..g7
  std::array<scalarvalueType,8> gvals;
  std::array<scalargradType,8>  ggrads;

  for (int i_g = 0; i_g < 8; i_g++)
    {
      gvals[i_g]  = variable_list.get_scalar_value(8 + i_g);
      ggrads[i_g] = variable_list.get_scalar_gradient(8 + i_g);
    }

  // Let's define placeholders for the grain constants (or read from your model constants).
  // e.g., double Wgrain=..., double Mgrain=..., double Kgrain=..., etc.
  // For demonstration, we re-use WV or define local placeholders:
  double Wgrain = WV;  // re-use your existing WV, or define new
  double Mgrain = 1.0; // or userInputs.get_model_constant_double("Mgrain"), etc.
  double Kgrain = 1.0e-14; // ...
  // no cross coupling for now

  for (int i_g = 0; i_g < 8; i_g++)
    {
      // polynomial derivative wrt g_i:  -W*g + W*g^3
      scalarvalueType Fp_gi = -Wgrain*gvals[i_g] + Wgrain*gvals[i_g]*gvals[i_g]*gvals[i_g];

      // The explicit residual update:
      //   g_i^{n+1} = g_i^n - dt*Mgrain*Fp_gi
      // And the gradient term => - dt*Mgrain*Kgrain * grad(g_i)
      scalarvalueType r_gval = gvals[i_g] - delt*Mgrain*Fp_gi;
      scalargradType  r_ggrd = -(delt*Mgrain*Kgrain) * ggrads[i_g];

      // Submit them:
      variable_list.set_scalar_value_term_RHS(8 + i_g, r_gval);
      variable_list.set_scalar_gradient_term_RHS(8 + i_g, r_ggrd);
    }
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// We keep the original code. We do not add anything for g0..g7 because they are
// explicit time-dependent fields. No changes needed.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // Unmodified original corrosion logic for mu, mupsi, Phi, irxn, etc.
  // ...
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// Again, no changes needed for the new grain fields, as they are not time-independent.

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // Original code to handle potential (Phi) LHS, no changes for the grains
  // ...
}

// =================================================================================
// thresholdField: a function particular to this app
// =================================================================================
// Remains unchanged:
template <int dim, int degree>
void
customPDE<dim, degree>::capFields(VectorizedArray<double> &ncp,
                                  VectorizedArray<double> &psicp,
                                  VectorizedArray<double>  n,
                                  VectorizedArray<double>  psi) const
{
  // Original method, unchanged
  // ...
}


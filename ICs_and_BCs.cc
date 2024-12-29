// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================
template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  // ---------------------------------------------------------------------
  // Original code for corrosion variables (indices 0..7)
  // ---------------------------------------------------------------------
  double epssqV = userInputs.get_model_constant_double("epssqV");
  double deltaV = std::sqrt(2.0 * epssqV);
  double posx   = p[0];
  double posy   = p[1];
  double cx     = 0.5 * userInputs.domain_size[0];
  double cy     = userInputs.domain_size[1];
  double rad    = std::sqrt((posx - cx) * (posx - cx) + (posy - cy) * (posy - cy));
  double n0pro  = 0.5 * (1.0 - std::tanh((rad0 - rad) / deltaV));

  if (index == 0) // "n"
    {
      scalar_IC = n0pro;
      if (scalar_IC > 1.0)
        scalar_IC = 1.0;
    }
  else if (index == 2) // "psi"
    {
      scalar_IC = 1.0 - n0pro;
      if (scalar_IC > 1.0)
        scalar_IC = 1.0;
    }
  else if (index == 5) // "cP"
    {
      scalar_IC = 1000.0;
    }
  else if (index == 6) // "Phi"
    {
      scalar_IC = 0.0;
    }
  // If it's one of the other original variables (1->mu,3->mupsi,4->cM,7->irxn),
  // default them to zero:
  else if (index < 8)
    {
      scalar_IC = 0.0;
    }
  // ---------------------------------------------------------------------
  // New logic for the 8 grain fields (indices 8..15)
  // ---------------------------------------------------------------------
  // We'll place 8 circular "grains" in the domain. For index=8 => g0, index=9 => g1, etc.
  else if (index >= 8 && index < 16)
    {
      // i_g is the grain index in 0..7
      unsigned int i_g = index - 8;

      // Define centers (in [0..1] fractional coordinates) and radii (fraction of domain size):
      static const std::vector<Point<2>> grain_centers =
      {
        Point<2>(0.25, 0.25),  // g0
        Point<2>(0.75, 0.25),  // g1
        Point<2>(0.25, 0.70),  // g2
        Point<2>(0.75, 0.70),  // g3
        Point<2>(0.50, 0.50),  // g4
        Point<2>(0.20, 0.50),  // g5
        Point<2>(0.80, 0.50),  // g6
        Point<2>(0.50, 0.15)   // g7
      };
      static const std::vector<double> grain_radii =
      {
        0.15,  // g0
        0.15,  // g1
        0.12,  // g2
        0.12,  // g3
        0.10,  // g4
        0.08,  // g5
        0.08,  // g6
        0.06   // g7
      };

      // Convert fractional center & radius into actual domain dimensions:
      double gx = grain_centers[i_g][0] * userInputs.domain_size[0];
      double gy = grain_centers[i_g][1] * userInputs.domain_size[1];
      double r  = grain_radii[i_g]      * userInputs.domain_size[0]; 
        // using X-scale for radius

      // distance from (posx,posy) to (gx,gy)
      double dist = std::sqrt( (posx - gx)*(posx - gx) + (posy - gy)*(posy - gy) );

      // For a smooth interface, use a tanh ramp of thickness 0.5 (arbitrary)
      double thickness = 0.5;
      scalar_IC = 0.5 * (1.0 - std::tanh( (dist - r)/ thickness ));
    }
  // everything else -> 0
  else
    {
      scalar_IC = 0.0;
    }
  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================
template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(
  [[maybe_unused]] const Point<dim>  &p,
  [[maybe_unused]] const unsigned int index,
  [[maybe_unused]] const unsigned int direction,
  [[maybe_unused]] const double       time,
  [[maybe_unused]] double            &scalar_BC,
  [[maybe_unused]] Vector<double>    &vector_BC)
{
  // --------------------------------------------------------------------------
  // Original code: left blank or with an optional snippet
  // --------------------------------------------------------------------------
  /*
if (index == 2){
  if (direction == 1){
      double x=p[0];
      double y=p[1];
      scalar_BC=std::sin(y/7.0);
  }
}
  */
  // -------------------------------------------------------------------------
}


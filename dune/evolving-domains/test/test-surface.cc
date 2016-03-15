/**************************************************************************

  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2014 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger
  Copyright (C) 2004 - 2014 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2014 Mirko Kraenkel
  Copyright (C) 2006 - 2014 Christoph Gersbacher
  Copyright (C) 2006 - 2014 Martin Nolte
  Copyright (C) 2011 - 2014 Tobias Malkmus
  Copyright (C) 2012 - 2014 Stefan Girke
  Copyright (C) 2013 - 2014 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner


  The dune-fem module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

**************************************************************************/
#include <config.h>

// iostream includes
#include <iostream>

// dune grid includes
#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

// refinement includes
#include <dune/fem/space/common/adaptmanager.hh>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include eoc output
#include <dune/fem/misc/femeoc.hh>
#include <dune/evolving-domains/gridwidth.hh>

// include geometrty grid part
#if USE_OLD_GRIDPART
#include <dune/fem/gridpart/geogridpart.hh>
#else
#include <dune/evolving-domains/gridpart/geogridpart.hh>
#endif
// include description of surface deformation
#include "deformation.hh"

// quadrature for area calculation
#include <dune/fem/quadrature/cachingquadrature.hh>

// timer
#include <ctime>

template< class GridPart, int polorder = 1 >
double computeArea( const GridPart& gridPart )
{
  // compute area using quadrature rules
  double ret = 0;
  typedef Dune::Fem::CachingQuadrature< GridPart, 0 > QuadratureType;
  typedef typename GridPart::GridViewType GridViewType;
  for( const auto& e : elements( static_cast<GridViewType>( gridPart ) ) )
    {
      const auto &geometry = e.geometry();

      QuadratureType quadrature( e, polorder+1 );
      const size_t numQuadraturePoints = quadrature.nop();

      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  //! [Compute local contribution of operator]
	  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

	  ret += weight;
	}
    }

  return ret;
}

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType, int polorder = 1 >
void algorithm ( HGridType &grid, int step, const int eocId )
{
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 1 > FunctionSpaceType;

  // we want to solve the problem on the leaf elements of the grid
  //! [Setup the grid part for a deforming domain]
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > HostGridPartType;
  HostGridPartType hostGridPart( grid );
  typedef DeformationCoordFunction< HGridType::dimensionworld > DeformationType;
  DeformationType deformation;

  typedef DiscreteDeformationCoordHolder< DeformationType, HostGridPartType, polorder > DiscreteDeformationCoordHolderType;
  typedef typename DiscreteDeformationCoordHolderType :: DiscreteFunctionType CoordinateFunctionType;
  DiscreteDeformationCoordHolderType discreteDeformation( deformation, hostGridPart );

  typedef Dune::Fem::GeoGridPart< CoordinateFunctionType > GridPartType;
  GridPartType gridPart( discreteDeformation.coordFunction() );
  //! [Setup the grid part for a deforming domain]

  // start timer
  const std::clock_t start = std::clock();

  // compute area
  const double area = computeArea< GridPartType, polorder >( gridPart );

  // stop timer
  const double timeEllapsed = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

  // compute error
  const double error = std::abs( area - 4.0 * M_PI );
  std::cout << "area:\t" << area << std::endl;

  // get errors
  std::vector< double > store = { error };
  Dune :: Fem :: FemEoc :: setErrors( eocId, store );

  // write to file / output
  const double h = EvolvingDomain :: GridWidth :: gridWidth( gridPart );
  const int dofs = discreteDeformation.dofs();
  Dune::Fem::FemEoc::write( h, dofs, timeEllapsed, 0.0, std::cout );
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file
  Dune::Fem::Parameter::append( "parameter" );

  // initialzie eoc file
  std::string eocOutPath = Dune::Fem::Parameter::getValue<std::string>("fem.eocOutPath", std::string("."));
  Dune::Fem::FemEoc::initialize( eocOutPath, "eoc", "surface only" );

  // add entries to eoc calculation
  std::vector<std::string> femEocHeaders = { "area error" };

  // get eoc id
  const int eocId = Dune::Fem::FemEoc::addEntry( femEocHeaders );

  // type of hierarchical grid
  typedef Dune :: AlbertaGrid< 2, 3 > HGridType;
  static_assert( HGridType :: dimension == HGridType :: dimensionworld - 1, "this code is written with the assumption grid dim = world dim -1" );

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "areatest.repeats", 0 );

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "areatest.level" );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // calculate first step
  algorithm< HGridType, POLORDER >( grid, 0, eocId );

  for( int step = 0; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    algorithm< HGridType, POLORDER >( grid, step, eocId );
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}

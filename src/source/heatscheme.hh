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
#ifndef HEAT_FEMSCHEME_HH
#define HEAT_FEMSCHEME_HH

// lagrange interpolation 
#include <dune/fem/operator/lagrangeinterpolation.hh>
// time provider
#include <dune/fem/solver/timeprovider.hh>

// local includes
#include "femscheme.hh" 

template< class GridPart >
struct ErrorOutput
{
  ErrorOutput( const GridPart& gridPart,
	       const Dune::Fem::TimeProviderBase &tp,
	       const DataOutputParameters &parameter )
    : gridPart_( gridPart ),
      tp_( tp ),
      maxh_( 0.0 ),
      maxTau_( 0.0 )
  {
    if( Dune::Fem::MPIManager::rank() == 0 )
      init( parameter );
  }

  ~ErrorOutput()
  {
    if( file_ and Dune::Fem::MPIManager::rank() == 0 )
      {
        file_ << "# h: " << maxh_ << std::endl;
        file_ << "# tau: " << maxTau_ << std::endl;
	file_.close();
      }
  }

  void write( const double l2Error, const double h1Error )
  {
    if( file_ and Dune::Fem::MPIManager::rank() == 0 )
      file_ << tp_.time() << "  " << l2Error << "  " << h1Error << std::endl;

    if( Dune::Fem::MPIManager::rank() == 0 )
      computeMeshSize();
  }

protected:
  void init( const DataOutputParameters &parameter )
  {
    std::string name = Dune :: Fem ::Parameter :: commonOutputPath() + "/";
    // add prefix for data file
    name += parameter.prefix();
    name += ".txt";

    std::cout << "opening file: " << name << std::endl;
    file_.open( name.c_str() );
    if( !file_ )
      {
	std::cout << "could not write error file" << std::endl;
      }

    if( file_ )
      file_ << "# time  $L^2$ error  $H^1$ error" << std::endl;
  }

  void computeMeshSize()
  {
    for( auto e = gridPart_.template begin< 0 >();
         e != gridPart_.template end< 0 >(); ++e )
      {
        const auto geo = e.geometry();
        for( unsigned int i = 0; i < geo.corners(); ++i )
          {
            for( unsigned int j = 0; j < i; ++j )
              {
                const double dist = ( geo.corner( i ) - geo.corner( j ) ).two_norm();
                maxh_ = std::max( dist, maxh_ );
              }
          }
      }

    maxTau_ = std::max( tp_.deltaT(), maxTau_ );
  }

private:
  const GridPart& gridPart_;
  const Dune::Fem::TimeProviderBase &tp_;

  double maxh_;
  double maxTau_;

  mutable std::ofstream file_;
};


// HeatScheme 
//-----------

template < class ImplicitModel, class ExplicitModel > 
struct HeatScheme : public FemScheme<ImplicitModel>
{
  typedef FemScheme<ImplicitModel> BaseType;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::ModelType ImplicitModelType;
  typedef ExplicitModel ExplicitModelType;
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;

  typedef Dune::Fem::L2Norm< GridPartType > L2NormType;
  typedef Dune::Fem::H1Norm< GridPartType > H1NormType;

  typedef ErrorOutput< GridPartType > ErrorOutputType;

  HeatScheme( GridPartType &gridPart, 
              const ImplicitModelType& implicitModel,
              const ExplicitModelType& explicitModel,
	      const unsigned int step = 0 )
  : BaseType(gridPart, implicitModel),
    explicitModel_(explicitModel),
    explicitOperator_( explicitModel_, discreteSpace_ ),
    errorOutput_( gridPart, implicitModel_.timeProvider(), DataOutputParameters( step ) ),
    linftyl2Error_( 0 ),
    l2h1Error_( 0 )
  {}

  void prepare() 
  { 
    // start rhs timer
    Dune::FemTimer::start( rhsIdx_ );

    // apply constraints, e.g. Dirichlet contraints, to the solution 
    explicitOperator_.prepare( explicitModel_.dirichletBoundary(), solution_ );
    // apply explicit operator and also setup right hand side 
    explicitOperator_( solution_, rhs_ );
    // apply constraints, e.g. Dirichlet contraints, to the result 
    explicitOperator_.prepare( solution_, rhs_ );

    // stop rhs timer
    Dune::FemTimer::stop( rhsIdx_ );
  }

  void initialize () 
  {
     Dune::Fem::LagrangeInterpolation
          < typename ExplicitModelType::InitialFunctionType, DiscreteFunctionType > interpolation;
     interpolation( explicitModel_.initialFunction(), solution_ );
  }

  template< class GridExactSolution >
  void closeTimestep( const GridExactSolution &exact, const double deltaT )
  {
    // find l2 error
    L2NormType l2norm( gridPart_ );
    const double l2Error = l2norm.distance( exact, solution_ );
    linftyl2Error_ = std::max( linftyl2Error_, l2Error );

    // find h1 error
    H1NormType h1norm( gridPart_ );
    const double h1Error = h1norm.distance( exact, solution_ );
    l2h1Error_ = std::sqrt( l2h1Error_*l2h1Error_ + deltaT * h1Error * h1Error );

    // write to file
    errorOutput_.write( l2Error, h1Error );
  }

  double linftyl2Error() const
  {
    return linftyl2Error_;
  }

  double l2h1Error() const
  {
    return l2h1Error_;
  }

private:
  using BaseType::gridPart_;
  using BaseType::discreteSpace_;
  using BaseType::solution_;
  using BaseType::implicitModel_;
  using BaseType::rhs_;
  const ExplicitModelType &explicitModel_;
  typename BaseType::EllipticOperatorType explicitOperator_; // the operator for the rhs 

  ErrorOutputType errorOutput_;
  double linftyl2Error_;
  double l2h1Error_;

  using BaseType::rhsIdx_;
};

#endif // end #if HEAT_FEMSCHEME_HH

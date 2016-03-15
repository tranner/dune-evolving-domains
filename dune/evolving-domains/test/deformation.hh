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
#ifndef DEFORMATION_HH
#define DEFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/space/common/functionspace.hh>

// DeformationCoordFunction
// ------------------------

template< int dimWorld >
struct BoundaryProjection
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
    y /= y.two_norm();
  }
};

// include discrete function space
#include <dune/fem/space/lagrange.hh>
// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>
// lagrange interpolation 
#include <dune/fem/operator/lagrangeinterpolation.hh>

template< class Deformation, class GridPart, const unsigned int codim, const unsigned int polorder = 1 >
class DiscreteDeformationCoordHolder;

template< class Deformation, class GridPart, const unsigned int polorder >
class DiscreteDeformationCoordHolder< Deformation, GridPart, 0, polorder >
{
  typedef Deformation DeformationType;
  typedef GridPart GridPartType;

public:
  typedef typename DeformationType :: FunctionSpaceType FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
  // choose type of discrete function
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  DiscreteDeformationCoordHolder( Deformation& deformation, GridPart& gridPart )
    : deformation_( deformation ), gridPart_( gridPart ),
      discreteSpace_( gridPart ),
      coordFunction_( "deformation", discreteSpace_ )
  {
    interpolate();
  }

  const DiscreteFunctionType& coordFunction() const
  {
    return coordFunction_;
  }

  int dofs() const
  {
    return discreteSpace_.size();
  }

  int elements() const
  {
    int n = 0;
    for( auto it = discreteSpace_.begin(); it != discreteSpace_.end(); ++it )
      ++n;
    return n;
  }

protected:
  void interpolate()
  {
    typedef typename DiscreteFunctionType::DofType DofType;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    // set all DoFs to infinity
    const DofIteratorType dend = coordFunction_.dend();
    for( DofIteratorType dit = coordFunction_.dbegin(); dit != dend; ++dit )
      *dit = std::numeric_limits< DofType >::infinity();

    for( const auto& entity : discreteSpace_ )
      {
	const auto& lagrangePointSet
	  = discreteSpace_.lagrangePointSet( entity );
	const int nop = lagrangePointSet.nop();

	auto df_local = coordFunction_.localFunction( entity );

	// does element contain a boundary segment?
	const bool boundary = entity.hasBoundaryIntersections();

	// if not on boundary, map is identity
	if( not boundary )
	  {
	    const auto geometry = entity.geometry();

	    // assume point based local dofs
	    int k = 0;
	    for( int qp = 0; qp < nop; ++qp )
	      {
		// if the first DoF for this point is already valid, continue
		if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
		  {
		    const auto hatx = coordinate( lagrangePointSet[ qp ] );

		    // evaluate the function in the Lagrange point
		    typename DiscreteFunctionType::RangeType phi
		      = geometry.global( hatx );

		    // assign the appropriate values to the DoFs
		    for( int i = 0; i < dimRange; ++i, ++k )
		      df_local[ k ] = phi[ i ];
		  }
		else
		  k += dimRange;
	      }
	  }
	else
	  {
	    const auto geometry = entity.geometry();
	    std::vector< bool > isVertexOnBoundary( entity.subEntities( dimRange ) );
	    for( int i = 0; i < geometry.corners(); ++i )
	      {
		isVertexOnBoundary[i] = std::abs( geometry.corner(i).two_norm() - 1.0 ) < 1.0e-8;
	      }

	    // assume point based local dofs
	    const int nop = lagrangePointSet.nop();
	    int k = 0;
	    for( int qp = 0; qp < nop; ++qp )
	      {
		// if the first DoF for this point is already valid, continue
		if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
		  {
		    const auto hatx = coordinate( lagrangePointSet[ qp ] );

		    Dune::FieldVector< double, dimRange+1 > lambda;
		    lambda[0] = 1.0;
		    for( int i = 0; i < dimRange; ++i )
		      {
			lambda[i+1] = hatx[i];
			lambda[0] -= hatx[i];
		      }

		    double lambdaStar = 0.0;
		    for( int i = 0; i < dimRange+1; ++i )
		      {
			if( isVertexOnBoundary[i] )
			  {
			    lambdaStar += lambda[i];
			  }
		      }

		    typename DiscreteFunctionType::RangeType phi
		      = geometry.global( hatx );
		    if( lambdaStar > 1.0e-8 )
		      {
			Dune::FieldVector< double, dimRange > y(0);
			for( int i = 0; i < dimRange+1; ++i )
			  {
			    if( isVertexOnBoundary[i] )
			      y.axpy(  lambda[i] / lambdaStar, geometry.corner(i) );
			  }

			Dune::FieldVector< double, dimRange > p;
			deformation_.evaluate( y, p );

			phi.axpy( std::pow( lambdaStar, df_local.order()+2 ), ( p - y ) );
		      }

		    // assign the appropriate values to the DoFs
		    for( int i = 0; i < dimRange; ++i, ++k )
		      df_local[ k ] = phi[ i ];
		  }
		else
		  k += dimRange;
	      }
	  }
      }
  }

private:
  Deformation& deformation_;
  GridPart& gridPart_;
  DiscreteFunctionSpaceType discreteSpace_;
  DiscreteFunctionType coordFunction_;
};

template< class Deformation, class GridPart, const unsigned int polorder >
class DiscreteDeformationCoordHolder< Deformation, GridPart, 1, polorder >
{
  typedef Deformation DeformationType;
  typedef GridPart GridPartType;

public:
  typedef typename DeformationType :: FunctionSpaceType FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
  // choose type of discrete function
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  DiscreteDeformationCoordHolder( Deformation& deformation, GridPart& gridPart )
    : deformation_( deformation ), gridPart_( gridPart ),
      discreteSpace_( gridPart ),
      coordFunction_( "deformation", discreteSpace_ )
  {
    interpolate();
  }

  const DiscreteFunctionType& coordFunction() const
  {
    return coordFunction_;
  }

  int dofs() const
  {
    return discreteSpace_.size();
  }

  int elements() const
  {
    int n = 0;
    for( auto it = discreteSpace_.begin(); it != discreteSpace_.end(); ++it )
      ++n;
    return n;
  }

protected:
  void interpolate()
  {
    Dune::Fem::LagrangeInterpolation
      < DeformationType, DiscreteFunctionType > interpolation;
    interpolation( deformation_, coordFunction_ );
  }

private:
  Deformation& deformation_;
  GridPart& gridPart_;
  DiscreteFunctionSpaceType discreteSpace_;
  DiscreteFunctionType coordFunction_;
};

#endif // #ifndef DEFORMATION_HH

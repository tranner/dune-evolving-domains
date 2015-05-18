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
struct DeformationCoordFunction
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;	

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  explicit DeformationCoordFunction ( const double time = 0.0 )
  : time_( time )
  {}

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    const double at = 2.0;// + 0.25 * sin( time_ );

    y[ 0 ] = x[ 0 ] * sqrt(at);
    y[ 1 ] = x[ 1 ];
    y[ 2 ] = x[ 2 ];
  }

  void setTime ( const double time ) { time_ = time; }

private:
  double time_;
};

//! deformation depending on a discrete function 
template <class DiscreteFunctionType>
class DeformationDiscreteFunction
: public Dune::DiscreteCoordFunction< double, 3, DeformationDiscreteFunction< DiscreteFunctionType > > 
{
  typedef Dune::DiscreteCoordFunction< double, 3, DeformationDiscreteFunction< DiscreteFunctionType > > BaseType;

  typedef typename DiscreteFunctionType :: GridType GridType ;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType ;
  typedef typename DiscreteFunctionType :: RangeType  RangeType ;
public:
  DeformationDiscreteFunction ( const DiscreteFunctionType& vertices )
  : vertices_( vertices )
  {}

  template< class HostEntity , class RangeVector >
  void evaluate ( const HostEntity &hostEntity, unsigned int corner,
                  RangeVector &y ) const
  { 
    DUNE_THROW(Dune::NotImplemented,"evaluate not implemented for codim > 0");
  }

  template <class RangeVector> 
  void evaluate ( const typename GridType :: template Codim<0>::Entity &entity, 
		  unsigned int corner, 
                  RangeVector &y ) const
  {  
    y = entity.geometry()[corner];

    return;
    typedef typename GridType::ctype  ctype;
    enum { dim = GridType::dimension };
   
    const Dune::ReferenceElement< ctype, dim > &refElement
      = Dune::ReferenceElements< ctype, dim >::general( entity.type() );

    LocalFunctionType localVertices = vertices_.localFunction( entity );	

    localVertices.evaluate( refElement.position( corner, dim ), y );
  }

  void setTime ( const double time )
  {
  }

protected:
  const DiscreteFunctionType& vertices_;
};

#endif // #ifndef DEFORMATION_HH


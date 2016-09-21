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
#ifndef POISSON_PROBLEMS_HH
#define POISSON_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "probleminterface.hh"

template <class FunctionSpace>
class SurfacePoissonProblem : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    phi[ 0 ] = 6.0 * x[0] * x[1];
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi[ 0 ] = x[0] * x[1];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret[ 0 ][ 0 ] = x[ 1 ] - 2.0 * x[ 0 ] * x[ 0 ] * x[ 1 ];
    ret[ 0 ][ 1 ] = x[ 0 ] - 2.0 * x[ 0 ] * x[ 1 ] * x[ 1 ];
    ret[ 0 ][ 2 ] = - 2.0 * x[ 0 ] * x[ 1 ] * x[ 2 ];
  }
  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
		 RangeType& value) const
  {
    value = RangeType( 0 );
  }
};

#endif // #ifndef POISSON_HH

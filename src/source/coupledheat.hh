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

#include "temporalprobleminterface.hh"

double Power( const double y, const int a )
{
  assert( a >= 0 );

  if( a == 0 )
    return y;

  return y * Power( y, a-1 );
}

template <class FunctionSpace>
class BulkHeatProblem
  : public TemporalProblemInterface< FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeFieldType       RangeFieldType;
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  BulkHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    phi = M_PI * cos( M_PI * time() ) * x[0] * x[1];
    phi *= alpha();
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0.0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to alpha * id
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = alpha();
  }

  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    // set to zero
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }

  //! robin coefficient
  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = alpha() * alpha() * RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( M_PI * time() ) * x[0] * x[1];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret = JacobianRangeType(0);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return false ;
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  virtual bool isNeumannPoint( const DomainType& x ) const
  {
    return true ;
  }

  // coupling data
  RangeFieldType alpha() const { return alpha_; }
  RangeFieldType beta() const { return beta_; }

private:
  const RangeFieldType alpha_;
  const RangeFieldType beta_;
};

template <class FunctionSpace>
class SurfaceHeatProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeFieldType       RangeFieldType;
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  SurfaceHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    phi = 3.0 * ( M_PI * cos( M_PI * time() ) + 6.0 * sin( M_PI * time() ) ) * x[0]*x[1];
    phi *= beta();
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0.0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to beta * id
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = beta();
  }

  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    // set to zero
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = beta() * beta() * RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = 3.0 * sin( M_PI * time() ) * x[0] * x[1];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret = JacobianRangeType(0);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return false ;
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  virtual bool isNeumannPoint( const DomainType& x ) const
  {
    return true ;
  }

  // coupling data
  RangeFieldType alpha() const { return alpha_; }
  RangeFieldType beta() const { return beta_; }

private:
  const RangeFieldType alpha_;
  const RangeFieldType beta_;
};

template <class FunctionSpace>
class ExchangeHeatProblem
  : public TemporalProblemInterface< FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeFieldType       RangeFieldType;
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  ExchangeHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {}

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to zero
    D = 0;
  }

  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    // set to zero
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    // set to zero
    m = -alpha() * beta() * RangeType(1);
  }

  virtual void d(const DomainType& x, RangeType &d) const
  {
    // set to zero
    d = 0;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return false ;
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  virtual bool isNeumannPoint( const DomainType& x ) const
  {
    return false ;
  }

  // coupling data
  RangeFieldType alpha() const { return alpha_; }
  RangeFieldType beta() const { return beta_; }

private:
  const RangeFieldType alpha_;
  const RangeFieldType beta_;
};


#endif // #ifndef POISSON_HH

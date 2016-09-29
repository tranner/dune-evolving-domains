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
  assert( a > 0 );

  if( a == 1 )
    return y;

  return y * Power( y, a-1 );
}

double Power( const double y, const double a )
{
  return std::pow( y, a );
}

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
    m = RangeType(0);
  }

  virtual void d(const DomainType& x, RangeType &d) const
  {
    // set to zero
    d = RangeType(0);
  }

  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = -alpha() * beta() * RangeType(1);
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

template <class FunctionSpace>
class NoExchangeHeatProblem
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

  NoExchangeHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider )
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
    m = RangeType(0);
  }

  virtual void d(const DomainType& x, RangeType &d) const
  {
    // set to zero
    d = RangeType(0);
  }

  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = RangeType(0);
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

template <class FunctionSpace, class DeformationType>
class BulkStationaryProblem
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

  BulkStationaryProblem( const Dune::Fem::TimeProviderBase &timeProvider,
			 const DeformationType& deformation )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) ),
      deformation_( deformation )
  {
    assert( std::abs( alpha_ - 1.0 ) < 1.0e-10 );
    assert( std::abs( beta_ - 1.0 ) < 1.0e-10 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    phi = cos( time() ) * x[1] * x[2];
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

  //! coefficient of time derivative
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = alpha() * RangeType(1);
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
    phi = sin( time() ) * x[1] * x[2];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret[0][0] = 0.0;
    ret[0][1] = 0.0;
    ret[0][2] = 0.0;
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

  const DeformationType& deformation_;
};

template <class FunctionSpace, class DeformationType >
class SurfaceStationaryProblem : public TemporalProblemInterface < FunctionSpace >
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

  SurfaceStationaryProblem( const Dune::Fem::TimeProviderBase &timeProvider,
			    const DeformationType& deformation )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) ),
      deformation_( deformation )
  {
    assert( std::abs( alpha_ - 1.0 ) < 1.0e-10 );
    assert( std::abs( beta_ - 1.0 ) < 1.0e-10 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double A = deformation_.a();
    const double q = sqrt( x[0] * x[0] / (A*A) + x[1]*x[1] + x[2]*x[2] );

    phi = (x[1]*x[2]*(-(Power(-1 + A,3)*Power(x[0],6)*((2 + q)*cos( time() ) + 2*sin( time() ))) +
		Power(A,4)*(Power(A,2)*(2 + q)*cos( time() ) +
			    2*(1 + A*(1 + q) + Power(A,2)*(5 + 2*q))*sin( time() )) +
		Power(-1 + A,2)*A*Power(x[0],4)*
		(3*A*(2 + q)*cos( time() ) + 2*(-3 + A*(7 + 2*q))*sin( time() )) -
		(-1 + A)*Power(A,2)*Power(x[0],2)*
		(3*Power(A,2)*(2 + q)*cos( time() ) +
		 2*(-3 + A*(-2 + q) + Power(A,2)*(11 + 4*q))*sin( time() ))))/
      (q*Power(Power(A,2) - (-1 + A)*Power(x[0],2),3));
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = RangeType(0);
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

  //! coefficient of time derivative
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = beta() * RangeType(1);
  }

  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    const double A = deformation_.a();
    const double q = sqrt( x[0] * x[0] / (A*A) + x[1]*x[1] + x[2]*x[2] );

    phi = x[1] * x[2] * sin( time() ) * ( 1.0 + 2.0 / q );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret = 0.0;
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

  const DeformationType& deformation_;
};

template <class FunctionSpace, class DeformationType >
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

  BulkHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider,
		   const DeformationType& deformation,
		   const bool coupled = true )
    : BaseType( timeProvider ),
      deformation_( deformation ),
      coupled_( coupled ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {
    assert( std::abs( alpha_ - 1.0 ) < 1.0e-10 );
    assert( std::abs( beta_ - 1.0 ) < 1.0e-10 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double a = deformation_.a();
    const double ap = deformation_.ap();

    phi = x[1]*x[2] * ( cos( time() ) + sin( time() ) * ap / ( 2.0 * a ) );
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    const double a = deformation_.a();
    const double q = x[0]*x[0] / (a*a) + x[1]*x[1] + x[2]*x[2];

    value = 2.0 * x[1]*x[2] * sin( time() ) / sqrt(q);

    if( coupled_ )
      value += alpha() * ( sin( time() ) * x[1]*x[2] )
	- beta() * ( sin( time() ) * x[0] * x[1] );
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

  //! coefficient of time derivative
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = alpha() * RangeType(1);
  }

  //! robin coefficient
  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = RangeType(0);
    if( coupled_ )
      a += alpha() * alpha() * RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( time() ) * x[1] * x[2];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret[0][0] = 0.0;
    ret[0][1] = sin( time() ) * x[2];
    ret[0][2] = sin( time() ) * x[1];
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
  const DeformationType &deformation_;
  const bool coupled_;
  const RangeFieldType alpha_;
  const RangeFieldType beta_;
};

template <class FunctionSpace, class DeformationType>
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

  SurfaceHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider,
		      const DeformationType& deformation,
		      const bool coupled = true )
    : BaseType( timeProvider ),
      deformation_( deformation ),
      coupled_( coupled ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {
    assert( std::abs( alpha_ - 1.0 ) < 1.0e-10 );
    assert( std::abs( beta_ - 1.0 ) < 1.0e-10 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    // define evolution of surface
    const double at = 1.0 + 0.25 * sin( time() );
    const double apt = 0.25 * cos( time() );

    // calculated surface parameters
    const double divGammaV = 0.5 * at * apt * ( x[1]*x[1] + x[2]*x[2] ) / ( x[0]*x[0] + at*at * ( x[1]*x[1] + x[2]*x[2] ) );
    const double N1 = 1/at * x[0] / sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] );
    const double N2 = x[1] / sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] );
    const double H = ( 2.0 * x[0] * x[0] + at * ( 1 + at ) * ( x[1]*x[1] + x[2]*x[2] ) )
      / ( sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] ) * ( x[0]*x[0] + at*at * ( x[1]*x[1] + x[2]*x[2] ) ) );


    // calculate solution and derivatives
    const double ux = sin( time() ) * x[0] * x[1];
    const double mdux = ( cos( time() ) + 0.5 * sin( time() ) * apt / at ) * x[0] * x[1];
    const double mlapux = sin( time() ) * (  2.0 * N1 * N2 + H * ( x[1] * N1 + x[0] * N2 ) );

    const double ubulkx = sin( time() ) * x[1] * x[2];

    // construct solution
    phi = mdux + divGammaV * ux + mlapux;
    if( coupled_ )
      phi += beta() * ux - alpha() * ubulkx;
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = RangeType(0);
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
    m = RangeType(0);
    if( coupled_ )
      m += beta() * beta() * RangeType(1);
  }

  //! coefficient of time derivative
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = beta() * RangeType(1);
  }

  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( time() ) * x[0] * x[1];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    // helper variables
    const double s = sin( time() );
    const double at = deformation_.a();

    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = s * x[1];
    grad[ 0 ][ 1 ] = s * x[0];
    grad[ 0 ][ 2 ] = 0.0;

    DomainType nu;
    nu[ 0 ] = 2.0 * x[ 0 ] / at;
    nu[ 1 ] = 2.0 * x[ 1 ];
    nu[ 2 ] = 2.0 * x[ 2 ];
    nu /= nu.two_norm();

    double dot = 0;
    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	dot += nu[ i ] * grad[ 0 ][ i ];
      }

    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	ret[ 0 ][ i ] = grad[ 0 ][ i ] - dot * nu[ i ];
      }
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
  const DeformationType deformation_;
  const bool coupled_;
  const RangeFieldType alpha_;
  const RangeFieldType beta_;
};

#endif // #ifndef POISSON_HH

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
		   const DeformationType& deformation )
    : BaseType( timeProvider ),
      deformation_( deformation ),
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
		      const DeformationType& deformation )
    : BaseType( timeProvider ),
      deformation_( deformation ),
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
    const double q = x[0]*x[0] / (a*a) + x[1]*x[1] + x[2]*x[2];
    const double t = time();

    auto Cos = []( double t ) -> double { return cos( t ); };
    auto Sin = []( double t ) -> double { return sin( t ); };
    auto Sqrt = []( double t ) -> double { return sqrt( t ); };

    phi = (x[1]*x[2]*(Power(a,7)*(4*Cos(t) + 20*Sin(t) + 2*Sqrt(q)*(Cos(t) + 4*Sin(t))) +
       2*Power(x[0],6)*Sin(t)*ap +
       a*Power(x[0],6)*(2*(2 + Sqrt(q))*Cos(t) + 4*Sin(t) -
	  (6 + Sqrt(q))*Sin(t)*ap) +
       Power(a,5)*(4*Sin(t) + Power(x[0],4)*
	   (6*(2 + Sqrt(q))*Cos(t) + 4*(7 + 2*Sqrt(q))*Sin(t)) +
	  Power(x[0],2)*(12*Cos(t) + 52*Sin(t) + 6*Sqrt(q)*(Cos(t) + 2*Sin(t)) -
	     3*(2 + Sqrt(q))*Sin(t)*ap)) +
       Power(a,3)*Power(x[0],2)*
	(-12*Sin(t) + Power(x[0],4)*
	   (6*((2 + Sqrt(q))*Cos(t) + 2*Sin(t)) -
	     (2 + Sqrt(q))*Sin(t)*ap) +
	  Power(x[0],2)*(12*Cos(t) + 52*Sin(t) + Sqrt(q)*(6*Cos(t) + 8*Sin(t)) -
	     4*(3 + Sqrt(q))*Sin(t)*ap)) +
       Power(a,2)*Power(x[0],4)*
	(Sin(t)*(-12 + (6 + Sqrt(q))*ap) +
	  Power(x[0],2)*(-6*((2 + Sqrt(q))*Cos(t) + 2*Sin(t)) +
	     2*(3 + Sqrt(q))*Sin(t)*ap)) +
       Power(a,4)*Power(x[0],2)*
	(-2*Power(x[0],4)*((2 + Sqrt(q))*Cos(t) + 2*Sin(t)) +
	  Power(x[0],2)*(-4*(6*Cos(t) + 17*Sin(t) + Sqrt(q)*(3*Cos(t) + 4*Sin(t))) +
	     3*(2 + Sqrt(q))*Sin(t)*ap) +
	  2*Sin(t)*(2 + 3*ap + Sqrt(q)*(2 + ap)))\
	+ Power(a,6)*(-2*Power(x[0],2)*
	   (6*Cos(t) + 22*Sin(t) + Sqrt(q)*(3*Cos(t) + 8*Sin(t))) +
	  Sin(t)*(4 + 2*ap + Sqrt(q)*(4 + ap)))))/
      (2.*Sqrt(q)*a*Power(Power(a,2) + Power(x[0],2) - a*Power(x[0],2),3));
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
    const double a = deformation_.a();
    const double q = x[0]*x[0] / (a*a) + x[1]*x[1] + x[2]*x[2];

    // return variable
    phi = sin( time() ) * x[1] * x[2] * ( 1 + 2 / sqrt(q) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    // helper variables
    const double s = sin( time() );
    const double f1 = std::sqrt( x[1]*x[1] + x[2]*x[2] + 16.0 * x[0]*x[0] / Power( 4.0 + s, 2 ) );
    const double fact = s * ( 1.0 + 2.0 / f1 );

    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = 0.0;
    grad[ 0 ][ 1 ] = fact * x[2];
    grad[ 0 ][ 2 ] = fact * x[1];

    const double at = 1.0 + 0.25 * sin( time() );
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
class BulkParabolicProblem
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

  BulkParabolicProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    // helper variables
    const double c = cos( time() );
    const double s = sin( time() );

    // return variable
    phi = 0.5 * x[1] * x[2] * c * ( 8 + 3 * s ) / ( 4 + s );
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
  const RangeFieldType alpha_;
  const RangeFieldType beta_;
};

template <class FunctionSpace>
class SurfaceParabolicProblem : public TemporalProblemInterface < FunctionSpace >
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

  SurfaceParabolicProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider ),
      alpha_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.alpha", 1.0 ) ),
      beta_( Dune::Fem::Parameter::getValue< RangeFieldType >( "coupled.beta", 1.0 ) )
  {
    assert(0);
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double tt = time();
    const double xx = x[0];
    const double yy = x[1];
    const double zz = x[2];

    const double x2 = xx*xx;
    const double y2 = yy*yy;
    const double z2 = zz*zz;

    const double c = cos(tt);
    const double s = sin(tt);

    const double g1 = yy*zz*c*(1 + (32*x2*s)/
			 (Power(4 + s,3)*Power(y2 + z2 + (16*x2)/Power(4 + s,2),1.5)) +
			 2/sqrt(y2 + z2 + (16*x2)/Power(4 + s,2)) -
			 (16*x2*s)/
			 (Power(4 + s,3)*Power((-4*x2*s + Power(4 + s,2))/Power(4 + s,2),1.5)));
    const double g2 = (yy*zz*(4 - 4*Power(xx,2) + sin(tt))*
     (2 + sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/
	 Power(4 + sin(tt),2)))*sin(2*tt))/
   (4.*Power(4 + sin(tt),2)*Power((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/
				  Power(4 + sin(tt),2),1.5));

    const double g3 = (yy*zz*(-5412528 - 650704*cos(2*tt) + 2688*cos(4*tt) -
       12288*Power(xx,4)*(177 + 2*cos(2*tt))*sin(tt) + 477088*Power(sin(tt),2) +
       648192*Power(xx,2)*Power(sin(tt),2) + 196608*Power(xx,4)*Power(sin(tt),2) -
       1024*Power(xx,4)*(905 + 6*cos(2*tt))*sin(tt)*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))\
	+ 247888*Power(sin(tt),2)*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))\
	+ 282368*Power(xx,2)*Power(sin(tt),2)*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))\
	+ 49152*Power(xx,4)*Power(sin(tt),2)*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))\
	- 46400*sin(3*tt) + 512*Power(xx,4)*
	(27860 + 1161*cos(2*tt) + 24*sin(3*tt)) -
       128*Power(xx,2)*sin(tt)*(22417 + 1257*cos(2*tt) + 24*sin(3*tt)) +
       128*Power(xx,4)*sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/
	  Power(4 + sin(tt),2))*(52059 + 1935*cos(2*tt) + 40*sin(3*tt)) -
       32*Power(xx,2)*sin(tt)*sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/
	  Power(4 + sin(tt),2))*(47557 + 2127*cos(2*tt) + 40*sin(3*tt)) +
       16*Power(xx,2)*(-550165 + 445*cos(2*tt) - 168*cos(4*tt) + 2132*sin(3*tt) -
	  4*sin(5*tt)) + 2*Power(yy,4)*(31 + cos(2*tt))*
	(88298 + 6440*cos(2*tt) - 42*cos(4*tt) + 725*sin(3*tt) - sin(5*tt)) +
       4*Power(yy,2)*Power(zz,2)*(31 + cos(2*tt))*
	(88298 + 6440*cos(2*tt) - 42*cos(4*tt) + 725*sin(3*tt) - sin(5*tt)) +
       2*Power(zz,4)*(31 + cos(2*tt))*
	(88298 + 6440*cos(2*tt) - 42*cos(4*tt) + 725*sin(3*tt) - sin(5*tt)) +
       Power(yy,4)*(31 + cos(2*tt))*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(92580 + 6832*cos(2*tt) - 44*cos(4*tt) + 773*sin(3*tt) - sin(5*tt)) +
       2*Power(yy,2)*Power(zz,2)*(31 + cos(2*tt))*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(92580 + 6832*cos(2*tt) - 44*cos(4*tt) + 773*sin(3*tt) - sin(5*tt)) +
       Power(zz,4)*(31 + cos(2*tt))*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(92580 + 6832*cos(2*tt) - 44*cos(4*tt) + 773*sin(3*tt) - sin(5*tt)) +
       64*sin(5*tt) + 8*sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/
	  Power(4 + sin(tt),2))*(-354827 - 42821*cos(2*tt) + 176*cos(4*tt) -
	  3092*sin(3*tt) + 4*sin(5*tt)) +
       4*sin(tt)*(138805 - 40669*cos(2*tt) + 168*cos(4*tt) - 2900*sin(3*tt) +
	  4*sin(5*tt)) - 8*Power(xx,2)*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(475045 - 8789*cos(2*tt) + 176*cos(4*tt) - 2452*sin(3*tt) + 4*sin(5*tt)) +
       8*Power(xx,2)*Power(yy,2)*sqrt((-4*Power(xx,2)*sin(tt) +
	    Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(1191611 + 69685*cos(2*tt) + 80*cos(4*tt) + 10*sin(tt) + 3712*sin(3*tt) +
	  6*sin(5*tt)) + 8*Power(xx,2)*Power(zz,2)*
	sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(1191611 + 69685*cos(2*tt) + 80*cos(4*tt) + 10*sin(tt) + 3712*sin(3*tt) +
	  6*sin(5*tt)) + sin(tt)*sqrt((-4*Power(xx,2)*sin(tt) +
	    Power(4 + sin(tt),2))/Power(4 + sin(tt),2))*
	(281898 - 85642*cos(2*tt) + 352*cos(4*tt) - 6184*sin(3*tt) + 8*sin(5*tt)) +
       16*Power(xx,2)*Power(yy,2)*
	(1233811 + 73501*cos(2*tt) + 144*cos(4*tt) + 12*sin(tt) + 3644*sin(3*tt) +
	  8*sin(5*tt)) + 16*Power(xx,2)*Power(zz,2)*
	(1233811 + 73501*cos(2*tt) + 144*cos(4*tt) + 12*sin(tt) + 3644*sin(3*tt) +
	  8*sin(5*tt))))/
   (8.*Power(-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2),3)*
    sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2)));

    const double g4 = (2*yy*zz*sin(tt))/
      sqrt((-4*Power(xx,2)*sin(tt) + Power(4 + sin(tt),2))/Power(4 + sin(tt),2));


    phi = g1+g2+g3+g4;
    assert( phi[0] == phi[0] );
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
    // helper variables
    const double s = sin( time() );
    const double f1 = std::sqrt( x[1]*x[1] + x[2]*x[2] + 16.0 * x[0]*x[0] / Power( 4.0 + s, 2 ) );

    // // return variable
    phi = x[1] * x[2] * s * ( 1.0 + 2.0 / f1 );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    // helper variables
    const double s = sin( time() );
    const double f1 = std::sqrt( x[1]*x[1] + x[2]*x[2] + 16.0 * x[0]*x[0] / Power( 4.0 + s, 2 ) );
    const double fact = s * ( 1.0 + 2.0 / f1 );

    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = 0.0;
    grad[ 0 ][ 1 ] = fact * x[2];
    grad[ 0 ][ 2 ] = fact * x[1];

    const double at = 1.0 + 0.25 * sin( time() );
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

#endif // #ifndef POISSON_HH

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

template <class FunctionSpace>
class SurfaceHeatProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
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
    : BaseType( timeProvider )
  {}

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

    // construct solution
    phi = mdux + divGammaV * ux + mlapux;
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0.0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to identity by default
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = 1;
  }

  //! advection coefficient (default = 0)
  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    // set to zero by default
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(1);
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
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * x[1];
    grad[ 0 ][ 1 ] = sin( time() ) * x[0];
    grad[ 0 ][ 2 ] = 0.0;

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
};

template <class FunctionSpace>
class SurfaceParabolicProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
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
    : BaseType( timeProvider )
  {}

  double Sin(const double a) const { return sin(a); }
  double Cos(const double a) const { return cos(a); }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double xx = x[0];
    const double yy = x[1];
    const double zz = x[2];
    const double tt = time();

    phi = (-12*Power(xx,4)*yy*Power(-1 + Cos(2*tt) - 8*Sin(tt),2) + 
     16*Power(xx,2)*yy*Sin(tt)*Power(4 + Sin(tt),4) - 
     yy*Sin(tt)*Power(4 + Sin(tt),5) + 
     4*Power(xx,5)*Power(Sin(tt),2)*
      (8*(48 + 16*Sin(tt) + Power(Sin(tt),2)) + 
        yy*(8*Cos(tt)*(3 + Sin(tt)) - 
           (4 + Sin(tt))*(128 + Power(Sin(tt),2) - 4*Sin(tt)*(-9 + Sin(xx*yy))))) + 
     (Power(xx,3)*Sin(tt)*Power(4 + Sin(tt),2)*
        (32*(-65 + Cos(2*tt) - 28*Sin(tt) + 2*Power(zz,2)*Sin(tt)*(6 + Sin(tt))) + 
          yy*(2656 - 96*Cos(2*tt) + Power(Sin(tt),3) - 32*Cos(tt)*(5 + 2*Sin(tt)) + 
             Sin(tt)*(1155 - 128*Sin(xx*yy)) - 16*Sin(xx*yy) - 
             16*Power(Sin(tt),2)*Sin(xx*yy) + 
             Power(Cos(tt),2)*(-3*Sin(tt) + 16*Sin(xx*yy)))))/4. + 
     (xx*Power(4 + Sin(tt),3)*(8*(-1 + Cos(2*tt) - 8*Sin(tt))*
           (-12 - Sin(tt) + 2*Power(zz,2)*(8 + Sin(tt))) + 
          yy*(32 + 132*Cos(tt) - 32*Cos(2*tt) - 4*Cos(3*tt) + 67*Cos(tt - xx*yy) - 
             Cos(3*tt - xx*yy) - 67*Cos(tt + xx*yy) + Cos(3*tt + xx*yy) + 
             506*Sin(tt) + 48*Sin(2*tt) + 2*Sin(3*tt) + 32*Sin(xx*yy) + 
             16*Sin(2*tt - xx*yy) - 16*Sin(2*tt + xx*yy))))/8.)/
      ((4 + Sin(tt))*Power(-4*Power(xx,2)*Sin(tt) + Power(4 + Sin(tt),2),2));

#if 0
    /* case b = 0 */
    phi = (xx*yy*(4*Power(xx,4)*Power(Sin(tt),2)*
        (8*Cos(tt)*(3 + Sin(tt)) - 
          (4 + Sin(tt))*(128 + Power(Sin(tt),2) - 4*Sin(tt)*(-9 + Sin(xx*yy)))) + 
       (Power(xx,2)*Sin(tt)*Power(4 + Sin(tt),2)*
          (2656 - 96*Cos(2*tt) + Power(Sin(tt),3) - 32*Cos(tt)*(5 + 2*Sin(tt)) + 
            Sin(tt)*(1155 - 128*Sin(xx*yy)) - 16*Sin(xx*yy) - 
            16*Power(Sin(tt),2)*Sin(xx*yy) + 
            Power(Cos(tt),2)*(-3*Sin(tt) + 16*Sin(xx*yy))))/4. + 
       Power(4 + Sin(tt),3)*(2*Cos(tt)*(8 + 6*Sin(tt) + Power(Sin(tt),2)) + 
          Sin(tt)*(Power(Sin(tt),2)*(-1 + Sin(xx*yy)) + 8*Sin(tt)*(1 + Sin(xx*yy)) + 
             16*(4 + Sin(xx*yy))))))/
      ((4 + Sin(tt))*Power(-4*Power(xx,2)*Sin(tt) + Power(4 + Sin(tt),2),2));
#endif
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0.0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to identity by default
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = 1;
    D *= 1 + x[0]*x[0];
  }

  //! advection coefficient (default = 0)
  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    AdvectionVectorType bt;
    bt[0] = 1;
    bt[1] = 2;
    bt[3] = 0;

    const double at = 1.0 + 0.25 * sin( time() );
    DomainType nu;
    nu[ 0 ] = 2.0 * x[ 0 ] / at;
    nu[ 1 ] = 2.0 * x[ 1 ];
    nu[ 2 ] = 2.0 * x[ 2 ];
    nu /= nu.two_norm();

    b = bt;
    nu *= ( bt * nu );
    b -= nu;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = sin( x[0]*x[1] );
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(1);
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
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * x[1];
    grad[ 0 ][ 1 ] = sin( time() ) * x[0];
    grad[ 0 ][ 2 ] = 0.0;

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
};

template <class FunctionSpace>
class BulkHeatProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
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
    : BaseType( timeProvider )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    // define evolution of bulk
    const double at = 1.0 + 0.25 * sin( time() );
    const double apt = 0.25 * cos( time() );
    const double t = time();

    phi = (sin(M_PI*x[1])*(2*at*(cos(t) + 2*Power(M_PI,2)*sin(t))*sin(M_PI*x[0]) +
			 sin(t)*(M_PI*cos(M_PI*x[0])*x[0] + sin(M_PI*x[0]))*apt))/(2.*at);
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    // define evolution of bulk
    const double at = 1.0 + 0.25 * sin( time() );
    const double t = time();

    value = (M_PI*sin(t)*(at*cos(M_PI*x[1])*x[1]*sin(M_PI*x[0]) + cos(M_PI*x[0])*x[0]*sin(M_PI*x[1])))/
      (at*sqrt(Power(x[0],2)/Power(at,2) + Power(x[1],2) + Power(x[2],2)));
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to identity by default
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = 1;
  }

  //! advection coefficient (default = 0)
  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    // set to zero by default
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( time() ) * sin( M_PI * x[0] ) * sin( M_PI * x[1] );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret[ 0 ][ 0 ] = sin( time() ) * M_PI * cos( M_PI * x[0] ) * sin( M_PI * x[1] );
    ret[ 0 ][ 1 ] = sin( time() ) * sin( M_PI * x[0] ) * M_PI * cos( M_PI * x[1] );
    ret[ 0 ][ 2 ] = 0.0;
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
};

template <class FunctionSpace>
class BulkParabolicProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
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
    : BaseType( timeProvider )
  {}

  double Sin(const double a) const { return sin(a); }
  double Cos(const double a) const { return cos(a); }
  double Sqrt( const double a ) const { return sqrt(a); }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double t = time();
    const double at = 1.0 + 0.25 * sin( time() );
    const double apt = 0.25 * cos( time() );

    phi = (4*Power(M_PI,2)*Cos(M_PI*x[0])*Cos(M_PI*x[1])*Power(x[0],2)*Sin(t) +
     (M_PI*Cos(M_PI*x[1])*x[0]*Sin(t)*Sin(M_PI*x[0])*(4*at - apt))/
      at + (2*at*(Cos(t)*Cos(M_PI*x[0])*Cos(M_PI*x[1]) +
           Sin(t)*(M_PI*Cos(M_PI*x[1])*Sin(M_PI*x[0]) +
              Cos(M_PI*x[0])*(Cos(M_PI*x[1])*(2*Power(M_PI,2) + Cos(x[0]*x[1])) +
                 2*M_PI*Sin(M_PI*x[1])))) +
	      Cos(M_PI*x[0])*Cos(M_PI*x[1])*Sin(t)*apt)/at)/2;
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    const double t = time();
    const double at = 1.0 + 0.25 * sin( time() );
    const double q = x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2];

    value = -((Sin(t)*(M_PI*Cos(M_PI*x[1])*Power(x[0],3)*Sin(M_PI*x[0]) + 
         Cos(M_PI*x[1])*x[0]*(-Cos(M_PI*x[0]) + M_PI*Sin(M_PI*x[0])) + 
         M_PI*at*Cos(M_PI*x[0])*Power(x[0],2)*x[1]*Sin(M_PI*x[1]) + 
		       at*Cos(M_PI*x[0])*x[1]*(-2*Cos(M_PI*x[1]) + M_PI*Sin(M_PI*x[1]))))/(Sqrt(q)*at));
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to identity by default
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = 1;
    D *= 1 + x[0]*x[0];
  }

  //! advection coefficient (default = 0)
  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    b[0] = 1;
    b[1] = 2;
    b[2] = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = cos( x[0]*x[1] );
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( time() ) * cos( M_PI * x[0] ) * cos( M_PI * x[1] );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    ret[ 0 ][ 0 ] = -sin( time() ) * M_PI * sin( M_PI * x[0] ) * cos( M_PI * x[1] );
    ret[ 0 ][ 1 ] = -sin( time() ) * cos( M_PI * x[0] ) * M_PI * sin( M_PI * x[1] );
    ret[ 0 ][ 2 ] = 0.0;
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
};

#endif // #ifndef POISSON_HH

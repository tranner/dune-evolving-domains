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

#include <dune/fem/io/parameter.hh>

#include "probleminterface.hh"

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
class BulkProblem
  : public ProblemInterface< FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeFieldType       RangeFieldType;
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  BulkProblem()
    : BaseType(),
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
    const double Xx = x[0], Yy = x[1];
    phi[ 0 ] =  beta()*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(3. + 4.*Xx - 4.*pow(Xx,2) + 4.*Yy - 4.*pow(Yy,2));
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
    m = alpha() * RangeType(1);
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
    phi[ 0 ] = beta() * exp( - x[0] * (x[0] - 1.0) - x[1] * (x[1] - 1.0)  );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    const double Xx = x[0], Yy = x[1];
    ret[0][0] = beta()*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. - 2.*Xx);
    ret[0][1] = beta()*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. - 2.*Yy);
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
};

template <class FunctionSpace>
class SurfaceProblem : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeFieldType       RangeFieldType;
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  SurfaceProblem()
    : BaseType(),
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
    const double Xx = x[0], Yy = x[1], Zz = x[2];
    phi[ 0 ] = beta()*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1.*Xx - 2.*pow(Xx,2) + 1.*Yy - 2.*pow(Yy,2)) + alpha()*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(3. + 4.*pow(Xx,5) - 4.*pow(Xx,6) + 4.*pow(Yy,5) - 4.*pow(Yy,6) + pow(Yy,4)*(13. - 4.*pow(Zz,2)) + pow(Xx,4)*(13. + 4.*Yy - 12.*pow(Yy,2) - 4.*pow(Zz,2)) + Yy*(8. - 2.*pow(Zz,2)) + pow(Yy,3)*(-10. + 4.*pow(Zz,2)) + pow(Xx,3)*(-10. - 2.*Yy + 8.*pow(Yy,2) + 4.*pow(Zz,2)) + pow(Yy,2)*(-14. + 5.*pow(Zz,2)) + pow(Xx,2)*(-14. + 8.*pow(Yy,3) - 12.*pow(Yy,4) + 5.*pow(Zz,2) + pow(Yy,2)*(26. - 8.*pow(Zz,2)) + Yy*(-10. + 4.*pow(Zz,2))) + Xx*(8. - 2.*pow(Yy,3) + 4.*pow(Yy,4) - 2.*pow(Zz,2) + Yy*(4. - 2.*pow(Zz,2)) + pow(Yy,2)*(-10. + 4.*pow(Zz,2)))) + pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(4. + 19.*Xx - 34.*pow(Xx,2) - 12.*pow(Xx,7) + 8.*pow(Xx,8) + 19.*Yy - 34.*pow(Yy,2) - 56.*pow(Yy,3) + 60.*pow(Yy,4) + 45.*pow(Yy,5) + 12.*Xx*pow(Yy,5) - 36.*pow(Xx,2)*pow(Yy,5) - 38.*pow(Yy,6) - 12.*Xx*pow(Yy,6) + 32.*pow(Xx,2)*pow(Yy,6) - 12.*pow(Yy,7) + 8.*pow(Yy,8) - 2.*Xx*pow(Zz,2) + 8.*pow(Xx,2)*pow(Zz,2) - 2.*Yy*pow(Zz,2) + 8.*pow(Yy,2)*pow(Zz,2) + 21.*pow(Yy,3)*pow(Zz,2) - 22.*pow(Yy,4)*pow(Zz,2) - 12.*pow(Yy,5)*pow(Zz,2) + 8.*pow(Yy,6)*pow(Zz,2) + pow(Xx,2)*pow(Yy,2)*(120. - 44.*pow(Zz,2)) + pow(Xx,2)*pow(Yy,3)*(88. - 24.*pow(Zz,2)) + Xx*pow(Yy,4)*(43. - 12.*pow(Zz,2)) + pow(Xx,5)*(45. + 12.*Yy - 36.*pow(Yy,2) - 12.*pow(Zz,2)) + Xx*Yy*(24. - 8.*pow(Zz,2)) + pow(Xx,6)*(-38. - 12.*Yy + 32.*pow(Yy,2) + 8.*pow(Zz,2)) + Xx*pow(Yy,3)*(-32. + 12.*pow(Zz,2)) + pow(Xx,2)*Yy*(-52. + 19.*pow(Zz,2)) + Xx*pow(Yy,2)*(-52. + 19.*pow(Zz,2)) + pow(Xx,2)*pow(Yy,4)*(-114. + 24.*pow(Zz,2)) + pow(Xx,3)*(-56. + 24.*pow(Yy,3) - 36.*pow(Yy,4) + 21.*pow(Zz,2) + pow(Yy,2)*(88. - 24.*pow(Zz,2)) + Yy*(-32. + 12.*pow(Zz,2))) + pow(Xx,4)*(60. - 36.*pow(Yy,3) + 48.*pow(Yy,4) - 22.*pow(Zz,2) + Yy*(43. - 12.*pow(Zz,2)) + pow(Yy,2)*(-114. + 24.*pow(Zz,2))));
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
    m = ( beta() * beta() + beta() ) * RangeType(1);
  }

  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = RangeType(0);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    const double Xx = x[0], Yy = x[1];
    phi[ 0 ] = exp(-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1.*alpha() + 1.*Xx - 2.*Xx*Xx + (1. - 2.*Yy)*Yy);
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    const double Xx = x[0], Yy = x[1];
    JacobianRangeType grad(0);
    grad[0][0] =  pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. + alpha()*(1. - 2.*Xx) - 4.*pow(Xx,2) + 4.*pow(Xx,3) + 1.*Yy - 2.*pow(Yy,2) + Xx*(-3. - 2.*Yy + 4.*pow(Yy,2)));
    grad[0][1] =  pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. + alpha()*(1. - 2.*Yy) + Xx*(1. - 2.*Yy) - 3.*Yy - 4.*pow(Yy,2) + 4.*pow(Yy,3) + pow(Xx,2)*(-2. + 4.*Yy));

    double dot = 0.0;
    for( int i = 0; i < dimDomain; ++i )
      dot += grad[ 0 ][ i ] * x[ i ];

    for( int i = 0; i < 3; ++i )
      ret[ 0 ][ i ] = grad[ 0 ][ i ] - dot * x[ i ];
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
class ExchangeProblem
  : public ProblemInterface< FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeFieldType       RangeFieldType;
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  ExchangeProblem()
    : BaseType(),
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

#endif // #ifndef POISSON_HH

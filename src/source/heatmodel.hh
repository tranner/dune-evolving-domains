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
#ifndef HEAT_MODEL_HH
#define HEAT_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

#include "temporalprobleminterface.hh"
#include "model.hh"

template< class FunctionSpace, class GridPart >
struct HeatModel : protected DiffusionModel<FunctionSpace,GridPart>
{
  typedef DiffusionModel<FunctionSpace,GridPart> BaseType;
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef TemporalProblemInterface< FunctionSpaceType > ProblemType ;

  typedef typename BaseType::ProblemType InitialFunctionType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  //! constructor taking problem reference, time provider, 
  //! time step factor( either theta or -(1-theta) ), 
  //! flag for the right hand side 
  HeatModel( const ProblemType& problem, 
             const GridPart &gridPart,
             const bool implicit ) 
    : BaseType(problem,gridPart),
      timeProvider_(problem.timeProvider()),
      implicit_( implicit ),
      timeStepFactor_( 0 )
  {
    // get theta for theta scheme 
    const double theta = Dune::Fem::Parameter::getValue< double >("heat.theta", 0.5 );
    if (implicit)
      timeStepFactor_ = theta ;
    else
      timeStepFactor_ = -( 1.0 - theta ) ;
  }

  template< class Entity, class Point >
  void source ( const Entity &entity, 
                const Point &x,
                const RangeType &value, 
                RangeType &flux ) const
  {
    linSource( value, entity, x, value, flux );
    // the explicit model should also evaluate the RHS
    if( !implicit_ ) 
    {
      const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
      // evaluate right hand side 
      RangeType rhs ;
      problem_.f( xGlobal, rhs );
      rhs  *= timeProvider_.deltaT();
      flux += rhs ;
    }
  }

  template< class Entity, class Point >
  void linSource ( const RangeType& uBar, 
                   const Entity &entity, 
                   const Point &x,
                   const RangeType &value, 
                   RangeType &flux ) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    RangeType m;
    problem_.m(xGlobal,m);
    for (unsigned int i=0;i<flux.size();++i)
      flux[i] = m[i]*value[i];
    flux *= timeStepFactor_ * timeProvider_.deltaT();
    // add term from time derivative
    flux += value;
  }

   //! return the diffusive flux 
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity, 
                       const Point &x,
                       const RangeType &value, 
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }
  //! return the diffusive flux 
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar, 
                          const JacobianRangeType &gradientBar,
                          const Entity &entity, 
                          const Point &x,
                          const RangeType &value, 
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
    flux  = gradient;
    flux *= timeStepFactor_ * timeProvider_.deltaT();
  }
  //! exact some methods from the problem class
  bool hasDirichletBoundary () const 
  {
    return BaseType::hasDirichletBoundary() ;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  bool isDirichletPoint( const DomainType& x ) const 
  {
    return BaseType::isDirichletPoint(x) ;
  }

  template< class Entity, class Point >
  void g( const RangeType& uBar, 
          const Entity &entity, 
          const Point &x,
          RangeType &u ) const
  {
    BaseType::g(uBar,entity,x,u);
  }

  // return Fem::Function for Dirichlet boundary values 
  typename BaseType::DirichletBoundaryType dirichletBoundary( ) const 
  {
    return BaseType::dirichletBoundary();
  }

  //! return reference to Problem's time provider 
  const TimeProviderType & timeProvider() const 
  {
    return timeProvider_;
  }

  const InitialFunctionType &initialFunction() const
  {
    return problem_;
  }

protected:
  using BaseType::problem_;
  const TimeProviderType &timeProvider_; 
  bool implicit_;
  double timeStepFactor_;

};
#endif // #ifndef HEAT_MODEL_HH

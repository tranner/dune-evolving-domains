// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef COUPLED_HEATSCHEME_HH
#define COUPLED_HEATSCHEME_HH

// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>
// time provider
#include <dune/fem/solver/timeprovider.hh>
// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// local includes
#include "coupledscheme.hh"

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
    init( parameter );
  }

  ~ErrorOutput()
  {
    if( file_ )
      {
        file_ << "# h: " << maxh_ << std::endl;
        file_ << "# tau: " << maxTau_ << std::endl;
      }

    if( file_ )
      file_.close();
  }

  void write( const double l2BulkError, const double h1BulkError,
              const double l2SurfaceError, const double h1SurfaceError )
  {
    if( file_ )
      file_ << tp_.time() << "  " << l2BulkError << "  " << h1BulkError
	    << "  " << l2SurfaceError << "  " << h1SurfaceError << std::endl;

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
      {
        file_ << "# time  $L^2(\\Omega(t))$ error  $H^1(\\Omega(t))$ error  "
              << "$L^2(\\Gamma(t))$ error  $H^1(\\Gamma(t))$ error" << std::endl;
      }
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

template< class ImplicitModel, class ExplicitModel, unsigned int codim >
struct TemporalFemSchemeHolder : public FemSchemeHolder< ImplicitModel, codim >
{
  typedef FemSchemeHolder< ImplicitModel, codim > BaseType;

  typedef typename BaseType :: GridType GridType;
  typedef typename BaseType :: GridPartType GridPartType;
  typedef typename BaseType :: ModelType ImplicitModelType;
  typedef ExplicitModel ExplicitModelType;
  typedef typename BaseType :: DiscreteFunctionType DiscreteFunctionType;

  TemporalFemSchemeHolder( GridPartType &gridPart,
                           const ImplicitModelType &implicitModel,
                           const ExplicitModelType &explicitModel )
    : BaseType( gridPart, implicitModel ),
      explicitModel_( explicitModel ),
      explicitOperator_( explicitModel, discreteSpace_ )
  {}

  TemporalFemSchemeHolder( const TemporalFemSchemeHolder &other )
    : BaseType( other ),
      explicitModel_( other.explicitModel_ ),
      explicitOperator_( other.explicitOperator_ )
  {}

  virtual void prepare()
  {
    // apply constraints, e.g. Dirichlet constraints, to the solution
    explicitOperator_.prepare( explicitModel_.dirichletBoundary(), solution() );
    // apply explicit operator and also setup right hand side
    explicitOperator_( solution(), rhs() );
    // apply constraints, e.g. Dirichlet constraints, to the result
    explicitOperator_.prepare( solution(), rhs() );
  }

  virtual void initialize ()
  {
     Dune::Fem::LagrangeInterpolation
          < typename ExplicitModelType::InitialFunctionType, DiscreteFunctionType > interpolation;
     interpolation( explicitModel_.initialFunction(), solution() );
  }

  using BaseType::gridPart;

private:
  using BaseType::gridPart_;
  using BaseType::discreteSpace_;
  using BaseType::implicitModel_;
  using BaseType::gridView;

  using BaseType::solution;
  using BaseType::rhs;

  const ExplicitModelType &explicitModel_;
  typename BaseType::EllipticOperatorType explicitOperator_; // the operator for the rhs
};

template< class BulkImplicitModel, class BulkExplicitModel,
	  class SurfaceImplicitModel, class SurfaceExplicitModel,
	  class ExchangeModel, class CoupledGrid >
class CoupledHeatScheme
  : public CoupledScheme< BulkImplicitModel, SurfaceImplicitModel,
			  ExchangeModel, CoupledGrid >
{
  typedef CoupledScheme< BulkImplicitModel, SurfaceImplicitModel,
			 ExchangeModel, CoupledGrid > BaseType;

public:
  typedef BulkImplicitModel BulkImplicitModelType;
  typedef BulkExplicitModel BulkExplicitModelType;
  typedef SurfaceImplicitModel SurfaceImplicitModelType;
  typedef SurfaceExplicitModel SurfaceExplicitModelType;
  typedef ExchangeModel ExchangeModelType;
  typedef CoupledGrid CoupledGridType;

  typedef typename BulkImplicitModelType :: GridPartType BulkGridPartType;
  typedef typename SurfaceImplicitModelType :: GridPartType SurfaceGridPartType;

  typedef TemporalFemSchemeHolder< BulkImplicitModelType, BulkExplicitModelType, 0 > BulkFemSchemeHolderType;
  typedef TemporalFemSchemeHolder< SurfaceImplicitModelType, SurfaceExplicitModelType, 1 > SurfaceFemSchemeHolderType;

  typedef typename BulkFemSchemeHolderType :: DiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename SurfaceFemSchemeHolderType :: DiscreteFunctionType SurfaceDiscreteFunctionType;

  typedef ErrorOutput< typename BulkFemSchemeHolderType :: GridPartType > ErrorOutputType;

  CoupledHeatScheme( BulkFemSchemeHolderType &bulkScheme, SurfaceFemSchemeHolderType &surfaceScheme,
                     ExchangeModelType &exchangeModel, const CoupledGridType &coupledGrid,
                     const unsigned int step = 0 )
    : BaseType( bulkScheme, surfaceScheme, exchangeModel, coupledGrid ),
      // error storage
      errorOutput_( bulkScheme.gridPart(), exchangeModel.timeProvider(),
                    DataOutputParameters( step ) ),
      linftyl2BulkError_( 0 ),
      l2h1BulkError_( 0 ),
      linftyl2SurfaceError_( 0 ),
      l2h1SurfaceError_( 0 )
  {}

  void prepare()
  {
    // start rhs timer
    Dune::FemTimer::start( rhsIdx_ );

    bulk().prepare();
    surface().prepare();

    // stop rhs timer
    Dune::FemTimer::stop( rhsIdx_ );
  }

  void initialize()
  {
    bulk().initialize();
    surface().initialize();
  }

  template< class BulkGridExactSolution, class SurfaceGridExactSolution,
            class TimeProvider >
  void closeTimestep( const BulkGridExactSolution &bulkExact,
                      const SurfaceGridExactSolution &surfaceExact,
                      const TimeProvider &tp, const bool first = false )
  {
    const double deltaT = tp.deltaT();

    const double bulkl2error = bulk().l2Error( bulkExact );
    linftyl2BulkError_ = std::max( linftyl2BulkError_, bulkl2error );

    const double bulkh1error = bulk().h1Error( bulkExact );
    l2h1BulkError_ = std::sqrt( l2h1BulkError_ * l2h1BulkError_ + deltaT * bulkh1error * bulkh1error );

    const double surfacel2error = surface().l2Error( surfaceExact );
    linftyl2SurfaceError_ = std::max( linftyl2SurfaceError_, surfacel2error );

    const double surfaceh1error = surface().h1Error( surfaceExact );
    l2h1SurfaceError_ = std::sqrt( l2h1SurfaceError_ * l2h1SurfaceError_ + deltaT * surfaceh1error * surfaceh1error );

    // write to file
    errorOutput_.write( bulkl2error, bulkh1error, surfacel2error, surfaceh1error );

#if 0
    // compute mass
    static double oldBulkMass = 0;
    const double bulkMass = bulk().mass();

    static double oldSurfaceMass = 0;
    const double surfaceMass = surface().mass();

    if( not first )
      {
        const double change = ( bulkMass - oldBulkMass )
          + ( surfaceMass - oldSurfaceMass );
        if( std::abs( change ) > std::sqrt(solverEps()) )
          {
            std::cout << "bulk mass: " << bulkMass << " change: " << ( bulkMass - oldBulkMass ) << "\n"
                      << "surface mass: " << surfaceMass << " change: " << ( surfaceMass - oldSurfaceMass ) << std::endl;
            std::cout << "sum change: " << change << std::endl;
            assert(0);
          }
      }

    oldBulkMass = bulkMass;
    oldSurfaceMass = surfaceMass;
#endif
  }

  double linftyl2BulkError() const
  {
    return linftyl2BulkError_;
  }
  double l2h1BulkError() const
  {
    return l2h1BulkError_;
  }
  double linftyl2SurfaceError() const
  {
    return linftyl2SurfaceError_;
  }
  double l2h1SurfaceError() const
  {
    return l2h1SurfaceError_;
  }

  using BaseType::bulkSolution;
  using BaseType::surfaceSolution;

protected:
  using BaseType::bulk;
  using BaseType::surface;

private:
  ErrorOutputType errorOutput_;
  double linftyl2BulkError_;
  double l2h1BulkError_;
  double linftyl2SurfaceError_;
  double l2h1SurfaceError_;

  using BaseType::solverEps;
  using BaseType::rhsIdx_;
};

#endif // #ifndef COUPLED_HEATSCHEME_HH

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

struct ErrorOutput
{
  ErrorOutput( const Dune::Fem::TimeProviderBase &tp,
	       const DataOutputParameters &parameter )
    : tp_( tp )
  {
    init( parameter );
  }

  ~ErrorOutput()
  {
    if( file_ )
      file_.close();
  }

  void write( const double l2BulkError, const double h1BulkError,
	      const double l2SurfaceError, const double h1SurfaceError )
  {
    if( file_ )
      file_ << tp_.time() << "  " << l2BulkError << "  " << h1BulkError
	    << "  " << l2SurfaceError << "  " << h1SurfaceError << std::endl;
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
      file_ << "# time  $L^2(\\Omega(t))$ error  $H^1(\\Omega(t))$ error  "
	    << "$L^2(\\Gamma(t))$ error  $H^1(\\Gamma(t))$ error" << std::endl;
  }

private:
  const Dune::Fem::TimeProviderBase &tp_;
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

  typedef Dune::Fem::L2Norm< GridPartType > L2NormType;
  typedef Dune::Fem::H1Norm< GridPartType > H1NormType;

  TemporalFemSchemeHolder( GridPartType &gridPart,
			   const ImplicitModelType &implicitModel,
			   const ExplicitModelType &explicitModel )
    : BaseType( gridPart, implicitModel ),
      explicitModel_( explicitModel ),
      explicitOperator_( explicitModel_ ),
      oldSolution_( solution_ )
  {}

  void prepare( bool firstCall )
  {
    if( firstCall )
      {
	oldSolution_.assign( solution_ );
      }

    // apply constraints, e.g. Dirichlet contraints, to the solution
    //explicitOperator_.prepare( explicitModel_.dirichletBoundary(), oldSolution_ );
    // apply explicit operator and also setup right hand side
    explicitOperator_( oldSolution_, rhs_ );
    // apply constraints, e.g. Dirichlet contraints, to the result
    //explicitOperator_.prepare( oldSolution_, rhs_ );
  }

  void initialize ()
  {
     Dune::Fem::LagrangeInterpolation
          < typename ExplicitModelType::InitialFunctionType, DiscreteFunctionType > interpolation;
     interpolation( explicitModel_.initialFunction(), solution_ );
  }

  template< class GridExactSolution >
  double l2Error( const GridExactSolution &exact ) const
  {
    L2NormType l2norm( gridPart_ );
    return l2norm.distance( exact, solution_ );
  }
  template< class GridExactSolution >
  double h1Error( const GridExactSolution &exact ) const
  {
    H1NormType h1norm( gridPart_ );
    return h1norm.distance( exact, solution_ );
  }

private:
  using BaseType::gridPart_;
  using BaseType::discreteSpace_;
  using BaseType::solution_;
  using BaseType::implicitModel_;
  using BaseType::rhs_;
  const ExplicitModelType &explicitModel_;
  typename BaseType::EllipticOperatorType explicitOperator_; // the operator for the rhs

  DiscreteFunctionType oldSolution_;
};

template< class CoupledImplicitModel, class CoupledExplicitModel,
	  class CoupledGrid >
class CoupledHeatScheme
{
public:
  typedef CoupledImplicitModel CoupledImplicitModelType;
  typedef CoupledExplicitModel CoupledExplicitModelType;

  typedef typename CoupledImplicitModelType :: BulkModelType BulkImplicitModelType;
  typedef typename CoupledImplicitModelType :: SurfaceModelType SurfaceImplicitModelType;

  typedef typename CoupledExplicitModelType :: BulkModelType BulkExplicitModelType;
  typedef typename CoupledExplicitModelType :: SurfaceModelType SurfaceExplicitModelType;

  typedef CoupledGrid CoupledGridType;

  typedef CoupledScheme< BulkImplicitModelType, SurfaceImplicitModelType, CoupledGridType >
  BaseType;

  typedef typename BulkImplicitModelType :: GridPartType BulkGridPartType;
  typedef typename SurfaceImplicitModelType :: GridPartType SurfaceGridPartType;

  typedef TemporalFemSchemeHolder< BulkImplicitModelType, BulkExplicitModelType, 0 > BulkFemSchemeHolderType;
  typedef TemporalFemSchemeHolder< SurfaceImplicitModelType, SurfaceExplicitModelType, 1 > SurfaceFemSchemeHolderType;

  typedef typename BulkFemSchemeHolderType :: DiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename SurfaceFemSchemeHolderType :: DiscreteFunctionType SurfaceDiscreteFunctionType;

  CoupledHeatScheme( BulkGridPartType &bulkGridPart,
		     SurfaceGridPartType &surfaceGridPart,
		     const CoupledImplicitModelType &coupledImplicitModel,
		     const CoupledExplicitModelType &coupledExplicitModel,
		     const CoupledGridType &coupledGrid,
		     const unsigned int step = 0 )
    : bulkScheme_( bulkGridPart, coupledImplicitModel.bulkModel(),
		   coupledExplicitModel.bulkModel() ),
      surfaceScheme_( surfaceGridPart, coupledImplicitModel.surfaceModel(),
		      coupledExplicitModel.surfaceModel() ),
    coupledImplicitModel_( coupledImplicitModel ),
    coupledExplicitModel_( coupledExplicitModel ),
    coupledGrid_( coupledGrid ),
    // tolerance for iterative solver
    solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) ),
    verbose_( Dune::Fem::Parameter::getValue< bool >( "coupled.solver.verbose", false ) ),
    // error stuffles
    errorOutput_( coupledImplicitModel.bulkModel().timeProvider(),
		  DataOutputParameters( step ) ),
    linftyl2BulkError_( 0 ),
    l2h1BulkError_( 0 ),
    linftyl2SurfaceError_( 0 ),
    l2h1SurfaceError_( 0 )
  {}

  const BulkDiscreteFunctionType &bulkSolution() const
  {
    return bulk().solution();
  }

  const SurfaceDiscreteFunctionType &surfaceSolution() const
  {
    return surface().solution();
  }
  void prepare( bool firstCall = true )
  {
    bulk().rhs().clear();
    surface().rhs().clear();

    bulk().prepare( firstCall );
    surface().prepare( firstCall );

    assembleRHS( coupledImplicitModel_,
		 bulkSolution(), surfaceSolution(),
		 coupledGrid_,
		 bulk().rhs(), surface().rhs() );
  }

  void initialize()
  {
    bulk().initialize();
    surface().initialize();
  }

  void solve( bool assemble )
  {
    //! [Solve the system]
    if( assemble )
    {
      // assemble linear operator (i.e. setup matrix)
      bulk().implicitOperator().jacobian( bulk().solution() , bulk().linearOperator() );
      surface().implicitOperator().jacobian( surface().solution() , surface().linearOperator() );
    }

    // inverse operator using linear operator
    typename BulkFemSchemeHolderType :: LinearInverseOperatorType bulkInvOp( bulk().linearOperator(), solverEps_, solverEps_ );
    typename SurfaceFemSchemeHolderType :: LinearInverseOperatorType surfaceInvOp( surface().linearOperator(), solverEps_, solverEps_ );

    const double eps = solverEps_ *
      ( bulk().rhs().scalarProductDofs( bulk().rhs() )
	+ surface().rhs().scalarProductDofs( surface().rhs() ) );

    iterations_ = 0;
    double update = 2.0 * eps;
    do
      {
	BulkDiscreteFunctionType oldBulkSolution( bulkSolution() );
	SurfaceDiscreteFunctionType oldSurfaceSolution( surfaceSolution() );

	// solve system
	bulkInvOp( bulk().rhs(), bulk().solution() );
	surfaceInvOp( surface().rhs(), surface().solution() );
	//! [Solve the system]

	// find difference
	oldBulkSolution -= bulkSolution();
	oldSurfaceSolution -= surfaceSolution();

	update = oldBulkSolution.scalarProductDofs( oldBulkSolution )
	  + oldSurfaceSolution.scalarProductDofs( oldSurfaceSolution );
	update = Dune::Fem::MPIManager::comm().sum( update );

	if( Dune::Fem::MPIManager::rank() == 0 and verbose_ )
	  std::cout << "it: " << iterations_ << " update: " << update
		    << " bulk solver it: " << bulkInvOp.iterations()
		    << " surface solver it: " << surfaceInvOp.iterations()
		    << std::endl;

	prepare( false );
	++iterations_;
      }
    while( update > eps );
  }

  template< class BulkGridExactSolution, class SurfaceGridExactSolution >
  void closeTimestep( const BulkGridExactSolution &bulkExact,
		      const SurfaceGridExactSolution &surfaceExact,
		      const double deltaT )
  {
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

  const int iterations() const
  {
    return iterations_;
  }

  const int dofs() const
  {
    return bulk().dofs() + surface().dofs();
  }

  const int elements() const
  {
    return bulk().elements() + surface().elements();
  }

protected:
  BulkFemSchemeHolderType &bulk() { return bulkScheme_; }
  const BulkFemSchemeHolderType &bulk() const { return bulkScheme_; }
  SurfaceFemSchemeHolderType &surface() { return surfaceScheme_; }
  const SurfaceFemSchemeHolderType &surface() const { return surfaceScheme_; }

private:
  BulkFemSchemeHolderType bulkScheme_;
  SurfaceFemSchemeHolderType surfaceScheme_;

  const CoupledImplicitModelType &coupledImplicitModel_;
  const CoupledExplicitModelType &coupledExplicitModel_;
  const CoupledGridType &coupledGrid_;

  const double solverEps_ ; // eps for linear solver
  const bool verbose_;
  unsigned int iterations_;

  ErrorOutput errorOutput_;
  double linftyl2BulkError_;
  double l2h1BulkError_;
  double linftyl2SurfaceError_;
  double l2h1SurfaceError_;
};

#endif // #ifndef COUPLED_HEATSCHEME_HH

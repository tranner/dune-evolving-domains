#ifndef COUPLED_HEATSCHEME_HH
#define COUPLED_HEATSCHEME_HH

// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>
// time provider
#include <dune/fem/solver/timeprovider.hh>

// local includes
#include "coupledscheme.hh"

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
      explicitOperator_( explicitModel_ )
  {}

  void prepare()
  {
    // apply constraints, e.g. Dirichlet contraints, to the solution
    explicitOperator_.prepare( explicitModel_.dirichletBoundary(), solution_ );
    // apply explicit operator and also setup right hand side
    explicitOperator_( solution_, rhs_ );
    // apply constraints, e.g. Dirichlet contraints, to the result
    explicitOperator_.prepare( solution_, rhs_ );
  }

  void initialize ()
  {
     Dune::Fem::LagrangeInterpolation
          < typename ExplicitModelType::InitialFunctionType, DiscreteFunctionType > interpolation;
     interpolation( explicitModel_.initialFunction(), solution_ );
  }

  template< class GridExactSolution >
  void closeTimestep( const GridExactSolution &exact, const double deltaT )
  {
  }

  double linftyl2Error() const
  {
    return 0.0;
  }

  double l2h1Error() const
  {
    return 0.0;
  }

private:
  using BaseType::gridPart_;
  using BaseType::discreteSpace_;
  using BaseType::solution_;
  using BaseType::implicitModel_;
  using BaseType::rhs_;
  const ExplicitModelType &explicitModel_;
  typename BaseType::EllipticOperatorType explicitOperator_; // the operator for the rhs
};

template< class BulkImplicitModel, class BulkExplicitModel,
	  class SurfaceImplicitModel, class SurfaceExplicitModel,
	  class CoupledGrid >
class CoupledHeatScheme
{
public:
  typedef BulkImplicitModel BulkImplicitModelType;
  typedef BulkExplicitModel BulkExplicitModelType;
  typedef SurfaceImplicitModel SurfaceImplicitModelType;
  typedef SurfaceExplicitModel SurfaceExplicitModelType;
  typedef CoupledGrid CoupledGridType;

  typedef CoupledScheme< BulkImplicitModel, SurfaceImplicitModel, CoupledGrid > BaseType;

  typedef typename BulkImplicitModelType :: GridPartType BulkGridPartType;
  typedef typename SurfaceImplicitModelType :: GridPartType SurfaceGridPartType;

  typedef TemporalFemSchemeHolder< BulkImplicitModelType, BulkExplicitModelType, 0 > BulkFemSchemeHolderType;
  typedef TemporalFemSchemeHolder< SurfaceImplicitModelType, SurfaceExplicitModelType, 1 > SurfaceFemSchemeHolderType;

  typedef typename BulkFemSchemeHolderType :: DiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename SurfaceFemSchemeHolderType :: DiscreteFunctionType SurfaceDiscreteFunctionType;

  CoupledHeatScheme( BulkGridPartType &bulkGridPart,
		     SurfaceGridPartType &surfaceGridPart,
		     const BulkImplicitModelType &bulkImplicitModel,
		     const BulkExplicitModelType &bulkExplicitModel,
		     const SurfaceImplicitModelType &surfaceImplicitModel,
		     const SurfaceExplicitModelType &surfaceExplicitModel,
		     const CoupledGridType &coupledGrid,
		     const unsigned int step = 0 )
    : bulkScheme_( bulkGridPart, bulkImplicitModel, bulkExplicitModel ),
      surfaceScheme_( surfaceGridPart, surfaceImplicitModel, surfaceExplicitModel ),
      coupledGrid_( coupledGrid ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) ),
      verbose_( Dune::Fem::Parameter::getValue< bool >( "coupled.solver.verbose", false ) )
  {}

  const BulkDiscreteFunctionType &bulkSolution() const
  {
    return bulk().solution();
  }

  const SurfaceDiscreteFunctionType &surfaceSolution() const
  {
    return surface().solution();
  }
  void prepare()
  {
    bulk().prepare();
    surface().prepare();

    assembleRHS( bulk().model(), surface().model(),
		 bulkSolution(), surfaceSolution(),
		 coupledGrid_,
		 bulk().rhs(), surface().rhs() );
#warning is forcing term in twice now?
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

	prepare();
	++iterations_;
      }
    while( update > eps );
  }

  template< class BulkGridExactSolution, class SurfaceGridExactSolution >
  void closeTimestep( const BulkGridExactSolution &bulkExact,
		      const SurfaceGridExactSolution &surfaceExplicitModel,
		      const double deltaT )
  {}

  double linftyl2BulkError() const
  {
    return 1.0;
  }
  double l2h1BulkError() const
  {
    return 1.0;
  }
  double linftyl2SurfaceError() const
  {
    return 1.0;
  }
  double l2h1SurfaceError() const
  {
    return 1.0;
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

  const CoupledGridType &coupledGrid_;

  const double solverEps_ ; // eps for linear solver
const bool verbose_;
  unsigned int iterations_;
};

#endif // #ifndef COUPLED_HEATSCHEME_HH

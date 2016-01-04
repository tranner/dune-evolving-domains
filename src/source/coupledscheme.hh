#ifndef COUPLED_FEMSCHEME_HH
#define COUPLED_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operator
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/linear/istloperator.hh>

// include solver
#include <dune/fem/solver/istlsolver.hh>

// local includes
#include "probleminterface.hh"
#include "rhs.hh"
#include "elliptic.hh"
#include "elliptic-coupled.hh"

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step, const std::string prefix = "" )
    : step_( step ), prefix_( prefix )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
    : step_( other.step_ ), prefix_( other.prefix_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "poisson-" + prefix_ + "-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
  const std::string prefix_;
};

template< class Model, unsigned int codim >
struct FemSchemeHolder
{
  //! type of the mathematical model
  typedef Model ModelType ;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of underlying grid view
  typedef typename GridPartType::GridViewType GridViewType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename ModelType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;

  //! define Laplace operator
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType > EllipticOperatorType;

  //! error help
  typedef Dune::Fem::L2Norm< GridPartType > L2NormType;
  typedef Dune::Fem::H1Norm< GridPartType > H1NormType;

  FemSchemeHolder( GridPartType &gridPart,
		   const ModelType& implicitModel)
    : implicitModel_( implicitModel ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      // the elliptic operator (implicit)
      implicitOperator_( implicitModel_, discreteSpace_ ),
      // create linear operator (domain space,range space)
      linearOperator_( "assembled elliptic operator", discreteSpace_, discreteSpace_ )
  {
    // set all DoF to zero
    solution_.clear();
  }

  FemSchemeHolder( const FemSchemeHolder& other )
    : implicitModel_( other.implicitModel_ ),
      gridPart_( other.gridPart_ ),
      discreteSpace_( other.gridPart_ ),
      solution_( other.solution_ ),
      rhs_( other.rhs_ ),
      implicitOperator_( other.implicitOperator_ ),
      linearOperator_( "assembled elliptic operator", discreteSpace_, discreteSpace_ )
  {}

  DiscreteFunctionType &solution() { return solution_; }
  const DiscreteFunctionType &solution() const { return solution_; }
  const ModelType &model()  const { return implicitModel_; }
  DiscreteFunctionType &rhs() { return rhs_; }
  const DiscreteFunctionType &rhs() const { return rhs_; }
  const DiscreteFunctionSpaceType &space() const { return discreteSpace_; }

  EllipticOperatorType &implicitOperator() { return implicitOperator_; }
  LinearOperatorType &linearOperator() { return linearOperator_; }

  const GridViewType gridView() const { return static_cast< GridViewType >( gridPart_ ); }

  //! set up the right hand side
  virtual void prepare()
  {
    // set boundary values for solution
    implicitOperator_.prepare( implicitModel_.dirichletBoundary(), solution_ );

    // assemble rhs
    assembleRHS ( implicitModel_.rightHandSide(), rhs_ );

    // apply constraints, e.g. Dirichlet contraints, to the result
    implicitOperator_.prepare( solution_, rhs_ );
  }

  virtual void initialize()
  {
    assert(0);
  }

  const int dofs() const
  {
    int tmp = discreteSpace_.size();
    return Dune::Fem::MPIManager::comm().sum( tmp );
  }

  const int elements() const
  {
    int elements_ = 0;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    const IteratorType end = discreteSpace_.end();
    for( IteratorType it = discreteSpace_.begin(); it != end; ++it )
      elements_++;
    return Dune::Fem::MPIManager::comm().sum( elements_ );
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

protected:
  const ModelType &implicitModel_; // the mathematical model

  GridPartType &gridPart_; // grid part(view), e.g. here the leaf grid the discrete space is build with

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_; // the unknown
  DiscreteFunctionType rhs_; // the right hand side

  EllipticOperatorType implicitOperator_; // the implicit operator

  LinearOperatorType linearOperator_; // the linear operator (i.e. jacobian of the implicit)
};

template< class BulkModel, class SurfaceModel, class ExchangeModel, class CoupledGrid >
class CoupledScheme
{
public:
  typedef BulkModel BulkModelType;
  typedef SurfaceModel SurfaceModelType;
  typedef ExchangeModel ExchangeModelType;
  typedef CoupledGrid CoupledGridType;

  typedef typename BulkModelType :: GridPartType BulkGridPartType;
  typedef typename SurfaceModelType :: GridPartType SurfaceGridPartType;

  typedef FemSchemeHolder< BulkModelType, 0 > BulkFemSchemeHolderType;
  typedef FemSchemeHolder< SurfaceModelType, 1 > SurfaceFemSchemeHolderType;

  typedef typename BulkFemSchemeHolderType :: DiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename SurfaceFemSchemeHolderType :: DiscreteFunctionType SurfaceDiscreteFunctionType;

  typedef Dune::Fem::SparseRowLinearOperator< BulkDiscreteFunctionType, SurfaceDiscreteFunctionType > BulkSurfaceLinearOperatorType;
  typedef Dune::Fem::SparseRowLinearOperator< SurfaceDiscreteFunctionType, BulkDiscreteFunctionType > SurfaceBulkLinearOperatorType;

  typedef DifferentiableMixingOperator< BulkSurfaceLinearOperatorType, ExchangeModelType, CoupledGridType >
    BulkSurfaceOperatorType;
  typedef DifferentiableMixingOperator< SurfaceBulkLinearOperatorType, ExchangeModelType, CoupledGridType >
    SurfaceBulkOperatorType;

#if 0
  CoupledScheme( BulkGridPartType &bulkGridPart, SurfaceGridPartType &surfaceGridPart,
		 const BulkModelType &bulkModel, const SurfaceModelType &surfaceModel,
		 const ExchangeModelType &exchangeModel, const CoupledGridType &coupledGrid )
    : bulkScheme_( bulkGridPart, bulkModel ),
      surfaceScheme_( surfaceGridPart, surfaceModel ),
      coupledGrid_( coupledGrid ),
      // the exchange operators
      bulkSurfaceImplicitOperator_( exchangeModel, bulk().space() ),
      surfaceBulkImplicitOperator_( exchangeModel, surface().space() ),
      // create exchange linear operators
      bulkSurfaceLinearOperator_( "assembled linear operator", bulk().space(), surface().space() ),
      surfaceBulkLinearOperator_( "assembled linear operator", surface().space(), bulk().space() ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) )
  {}
#endif

  CoupledScheme( BulkFemSchemeHolderType &bulkScheme, SurfaceFemSchemeHolderType &surfaceScheme,
		 const ExchangeModelType &exchangeModel, const CoupledGridType &coupledGrid )
  : bulkScheme_( bulkScheme ),
    surfaceScheme_( surfaceScheme ),
    coupledGrid_( coupledGrid ),
    // the exchange operators
    bulkSurfaceImplicitOperator_( exchangeModel, coupledGrid ),
    surfaceBulkImplicitOperator_( exchangeModel, coupledGrid ),
    // create exchange linear operators
    bulkSurfaceLinearOperator_( "assembled linear operator", bulk().space(), surface().space() ),
    surfaceBulkLinearOperator_( "assembled linear operator", surface().space(), bulk().space() ),
    // tolerance for iterative solver
    solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) ),
    maxIter_( Dune::Fem::Parameter::getValue< unsigned int >( "coupled.maxiter", 1000 ) ),
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
  }

  void solve( bool assemble )
  {
    //! [Solve the system]
    if( assemble )
    {
      // assemble linear operator (i.e. setup matrix)
      bulk().implicitOperator().jacobian( bulk().solution() , bulk().linearOperator() );
      surface().implicitOperator().jacobian( surface().solution() , surface().linearOperator() );

      bulkSurfaceImplicitOperator_.jacobian( bulk().solution(), bulkSurfaceLinearOperator_ );
      surfaceBulkImplicitOperator_.jacobian( surface().solution(), surfaceBulkLinearOperator_ );
    }

    // inverse operator using linear operator
    const double innerSolverEps = solverEps_ * 1.0e-2;
    typename BulkFemSchemeHolderType :: LinearInverseOperatorType bulkInvOp( bulk().linearOperator(), innerSolverEps, innerSolverEps );
    typename SurfaceFemSchemeHolderType :: LinearInverseOperatorType surfaceInvOp( surface().linearOperator(), innerSolverEps, innerSolverEps );

    const double eps = std::max( solverEps_ * solverEps_ *
      ( bulk().rhs().normSquaredDofs()+ surface().rhs().normSquaredDofs() ), 1.0e-16 );

    iterations_ = 0;
    double update = 2.0 * eps;
    for( iterations_ = 0; iterations_ < maxIter_; ++iterations_ )
      {
	// store old solution
	BulkDiscreteFunctionType oldBulkSolution( bulk().rhs() );
	SurfaceDiscreteFunctionType oldSurfaceSolution( surface().rhs() );

	// construct new rhs
	BulkDiscreteFunctionType myBulkRhs( bulk().rhs() );
	SurfaceDiscreteFunctionType mySurfaceRhs( surface().rhs() );

	BulkDiscreteFunctionType Bv( bulkSolution() );
	surfaceBulkLinearOperator_( surfaceSolution(), Bv );
	SurfaceDiscreteFunctionType Btu( surfaceSolution() );
	bulkSurfaceLinearOperator_( bulkSolution(), Btu );

	myBulkRhs -= Bv;
	mySurfaceRhs -= Btu;

	// solve system
	bulkInvOp( myBulkRhs, bulk().solution() );
	surfaceInvOp( mySurfaceRhs, surface().solution() );
	//! [Solve the system]

	// find difference
	oldBulkSolution -= bulkSolution();
	oldSurfaceSolution -= surfaceSolution();

	update = oldBulkSolution.scalarProductDofs( oldBulkSolution )
	  + oldSurfaceSolution.scalarProductDofs( oldSurfaceSolution );
	update = Dune::Fem::MPIManager::comm().sum( update );

	if( update < eps )
	  break;

	if( Dune::Fem::MPIManager::rank() == 0 and verbose_ )
	  {
	    std::cout << " it: " << iterations_
		    << " update: " << update << " / " << eps
		    << " bulk solver it: " << bulkInvOp.iterations()
		    << " surface solver it: " << surfaceInvOp.iterations()
		    << std::endl;
	    std::cout << oldBulkSolution.scalarProductDofs( oldBulkSolution ) << " and "
		      << oldSurfaceSolution.scalarProductDofs( oldSurfaceSolution ) << std::endl;

	    if( iterations_ > 5 )
	      { int qq; std::cin >> qq; }
	  }

	prepare();
      }
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

  const double solverEps() const { return solverEps_; }

private:
  BulkFemSchemeHolderType &bulkScheme_;
  SurfaceFemSchemeHolderType &surfaceScheme_;

  const CoupledGridType &coupledGrid_;

  const BulkSurfaceOperatorType bulkSurfaceImplicitOperator_;
  const SurfaceBulkOperatorType surfaceBulkImplicitOperator_;

  BulkSurfaceLinearOperatorType bulkSurfaceLinearOperator_;
  SurfaceBulkLinearOperatorType surfaceBulkLinearOperator_;

  const double solverEps_ ; // eps for linear solver
  unsigned int iterations_;
  const unsigned int maxIter_;
  const bool verbose_;
};

#endif // #ifndef COUPLED_FEMSCHEME_HH

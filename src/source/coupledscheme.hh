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

// include matrix
#include <dune/istl/bcrsmatrix.hh>

// include solver
#include <dune/fem/solver/istlsolver.hh>

// error norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// timer
#include <dune/fem/misc/femtimer.hh>

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
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
#if 0
  typedef Dune::Fem::ISTLGMResOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
#endif

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
    rhs_.clear();
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
  const LinearOperatorType &linearOperator() const { return linearOperator_; }

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

  const int nDofs() const
  {
    int tmp = discreteSpace_.size();
    return Dune::Fem::MPIManager::comm().sum( tmp );
  }

  const int nElements() const
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

  double mass() const
  {
    using QuadratureType = Dune::Fem::CachingQuadrature< GridPartType, 0 >;

    double ret = 0;

    // element loop
    for( auto e : elements( gridView() ) )
      {
        // find element information
        const auto geo = e.geometry();
        const auto& uLocal = solution().localFunction( e );

        // construct quadrature
        QuadratureType quadrature( e, discreteSpace_.order() );
        const unsigned int numQuadraturePoints = quadrature.nop();

        // quadrature loop
        for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
          {
            // obtain quadrature rule
            const auto& xLocal = quadrature.point( pt );
            const double weight = quadrature.weight( pt ) * geo.integrationElement( xLocal );

            // evaluate solution
            typename DiscreteFunctionType :: RangeType ux;
            uLocal.evaluate( quadrature[ pt ], ux );

            // multiply by quadrature weight
            ux *= weight;

            // add to result
            ret += ux;
          }
      }

    return ret;
  }

  const GridPartType& gridPart() const { return gridPart_; }

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

  // types of full system matrix
  using BlockType = typename BulkDiscreteFunctionType :: DofBlockType;
  using FullMatrixType = typename Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > >;

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
    verbose_( Dune::Fem::Parameter::getValue< bool >( "coupled.solver.verbose", false ) ),
    // timer indexes
    rhsIdx_( Dune::FemTimer::addTo( "rhs", 1 ) ),
    matrixIdx_( Dune::FemTimer::addTo( "matrix", 1 ) ),
    solverIdx_( Dune::FemTimer::addTo( "solver", 1 ) )
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
    // reset timers
    Dune::FemTimer::reset();

    // start rhs timer
    Dune::FemTimer::start( rhsIdx_ );

    bulk().prepare();
    surface().prepare();

    // stop rhs timer
    Dune::FemTimer::stop( rhsIdx_ );
  }

  void solve( bool assemble )
  {
    //! [Solve the system]
    if( assemble )
    {
      Dune::FemTimer::start( matrixIdx_ );

      // assemble linear operator (i.e. setup matrix)
      bulk().implicitOperator().jacobian( bulk().solution() , bulk().linearOperator() );
      surface().implicitOperator().jacobian( surface().solution() , surface().linearOperator() );

      bulkSurfaceImplicitOperator_.jacobian( bulk().solution(), bulkSurfaceLinearOperator_ );
      surfaceBulkImplicitOperator_.jacobian( surface().solution(), surfaceBulkLinearOperator_ );

      Dune::FemTimer::stop( matrixIdx_ );
    }

    Dune::FemTimer::start( solverIdx_ );

    // number of dofs for bulk and surface vectors
    const unsigned int nBulk = bulk().nDofs();
    const unsigned int nSurf = surface().nDofs();

    // construct matrix
    FullMatrixType systemMatrix( nBulk + nSurf, nBulk + nSurf, FullMatrixType :: random );
    assembleSystemMatrix( systemMatrix );

    using DF = Dune::BlockVector< BlockType >;

    DF rhs( nBulk + nSurf );
    for( unsigned int i = 0; i < nBulk; ++i )
      {
	rhs[ i ] = *bulk().rhs().block(i);
      }
    for( unsigned int i = 0; i < nSurf; ++i )
      {
	rhs[ i+nBulk ] = *surface().rhs().block(i);
      }

    Dune::MatrixAdapter< FullMatrixType, DF, DF > systemMatrixAdapter( systemMatrix );
    // Dune::SeqJac< FullMatrixType, DF, DF > preconditioner( systemMatrix, 1, 0.1 );
    Dune::SeqGS< FullMatrixType, DF, DF > preconditioner( systemMatrix, 1, 0.1 );
    Dune::SeqScalarProduct< DF > scp;

    using Solver = typename Dune::RestartedGMResSolver< DF >;
    Solver solver( systemMatrixAdapter, scp,
     		   preconditioner,
    		   solverEps_, 5,
		   std::numeric_limits< int >::max(),
    		   false );

    DF x( nBulk + nSurf );
    Dune::InverseOperatorResult res;
    solver.apply( x, rhs, res );

    for( unsigned int i = 0; i < nBulk; ++i )
      {
	*bulk().solution().block(i) = x[ i ];
      }
    for( unsigned int i = 0; i < nSurf; ++i )
      {
	*surface().solution().block(i) = x[ i+nBulk ];
      }

    Dune::FemTimer::stop( solverIdx_ );

    // extract results data
    iterations_ = res.iterations;
    if( not res.converged )
      std::cerr << "warning: solver not converged" << std::endl;
  }

  const int iterations() const
  {
    return iterations_;
  }

  const int nDofs() const
  {
    return bulk().nDofs() + surface().nDofs();
  }

  const int nElements() const
  {
    return bulk().nElements() + surface().nElements();
  }

  void printTimers() const
  {
    Dune::FemTimer::print( std::cout, "Timing data" );
    Dune::FemTimer::removeFrom( rhsIdx_ );
    Dune::FemTimer::removeFrom( matrixIdx_ );
    Dune::FemTimer::removeFrom( solverIdx_ );
  }

protected:
  void assembleSystemMatrix( FullMatrixType& systemMatrix ) const
  {
    // number of dofs for bulk and surface vectors
    const unsigned int nBulk = bulk().nDofs();
    const unsigned int nSurf = surface().nDofs();

    // extract matrices
    const auto& bulkBulkMatrix = bulk().linearOperator().systemMatrix().matrix();
    const auto& surfSurfMatrix = surface().linearOperator().systemMatrix().matrix();
    const auto& bulkSurfMatrix = surfaceBulkLinearOperator_.systemMatrix().matrix();
    const auto& surfBulkMatrix = bulkSurfaceLinearOperator_.systemMatrix().matrix();

    // initially set row size for each row
    for( unsigned int i = 0; i < nBulk; ++i )
      {
	const unsigned int bbNnzRow = bulkBulkMatrix.numNonZeros( i );
	const unsigned int bsNnzRow = bulkSurfMatrix.numNonZeros( i );
	systemMatrix.setrowsize( i, bbNnzRow + bsNnzRow );
      }
    for( unsigned int i = 0; i < nSurf; ++i )
      {
	const unsigned int sbNnzRow = surfBulkMatrix.numNonZeros( i );
	const unsigned int ssNnzRow = surfSurfMatrix.numNonZeros( i );
	systemMatrix.setrowsize( i + nBulk, sbNnzRow + ssNnzRow );
      }
    systemMatrix.endrowsizes();

    // add column entries to row
    for( unsigned int i = 0; i < nBulk; ++i )
      {
	const unsigned int bbNnzRow = bulkBulkMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < bbNnzRow; ++fakeCol )
	  {
	    const int j = bulkBulkMatrix.realCol( i, fakeCol );
	    systemMatrix.addindex( i, j );
	  }

	const unsigned int bsNnzRow = bulkSurfMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < bsNnzRow; ++fakeCol )
	  {
	    const int j = bulkSurfMatrix.realCol( i, fakeCol );
	    systemMatrix.addindex( i, j + nBulk );
	  }
      }
    for( unsigned int i = 0; i < nSurf; ++i )
      {
	const unsigned int sbNnzRow = surfBulkMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < sbNnzRow; ++fakeCol )
	  {
	    const int j = surfBulkMatrix.realCol( i, fakeCol );
	    systemMatrix.addindex( i + nBulk, j );
	  }

	const unsigned int ssNnzRow = surfSurfMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < ssNnzRow; ++fakeCol )
	  {
	    const int j = surfSurfMatrix.realCol( i, fakeCol );
	    systemMatrix.addindex( i + nBulk, j + nBulk );
	  }
      }

    // finalize column setup phase
    systemMatrix.endindices();

    // set entries using random access
    for( unsigned int i = 0; i < nBulk; ++i )
      {
	const unsigned int bbNnzRow = bulkBulkMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < bbNnzRow; ++fakeCol )
	  {
	    const auto p = bulkBulkMatrix.realValue( i, fakeCol );
	    const auto& val = p.first;
	    const unsigned int j = p.second;
	    systemMatrix[i][j] = val;
	  }

	const unsigned int bsNnzRow = bulkSurfMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < bsNnzRow; ++fakeCol )
	  {
	    const auto p = bulkSurfMatrix.realValue( i, fakeCol );
	    const auto& val = p.first;
	    const unsigned int j = p.second;
	    systemMatrix[i][j+nBulk] = val;
	  }
      }
    for( unsigned int i = 0; i < nSurf; ++i )
      {
	const unsigned int sbNnzRow = surfBulkMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < sbNnzRow; ++fakeCol )
	  {
	    const auto p = surfBulkMatrix.realValue( i, fakeCol );
	    const auto& val = p.first;
	    const unsigned int j = p.second;
	    systemMatrix[i+nBulk][j] = val;
	  }

	const unsigned int ssNnzRow = surfSurfMatrix.numNonZeros( i );
	for( unsigned int fakeCol = 0; fakeCol < ssNnzRow; ++fakeCol )
	  {
	    const auto p = surfSurfMatrix.realValue( i, fakeCol );
	    const auto& val = p.first;
	    const unsigned int j = p.second;
	    systemMatrix[i+nBulk][j+nBulk] = val;
	  }
      }
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

protected:
  unsigned int rhsIdx_;
  unsigned int matrixIdx_;
  unsigned int solverIdx_;
};

#endif // #ifndef COUPLED_FEMSCHEME_HH

#ifndef COUPLED_FEMSCHEME_HH
#define COUPLED_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operator
#include <dune/fem/operator/linear/istloperator.hh>

// include solver
#include <dune/fem/solver/istlsolver.hh>

// local includes
#include "probleminterface.hh"
#include "rhs-coupled.hh"
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

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename ModelType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;

  //! define Laplace operator
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType, codim > EllipticOperatorType;

  FemSchemeHolder( GridPartType &gridPart,
		   const ModelType& implicitModel)
    : implicitModel_( implicitModel ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      // the elliptic operator (implicit)
      implicitOperator_( implicitModel_ ),
      // create linear operator (domainSpace,rangeSpace)
      linearOperator_( "assembled elliptic operator", discreteSpace_, discreteSpace_ )
  {
    // set all DoF to zero
    solution_.clear();
  }

  DiscreteFunctionType &solution() { return solution_; }
  const DiscreteFunctionType &solution() const { return solution_; }
  const ModelType &model()  const { return implicitModel_; }
  DiscreteFunctionType &rhs() { return rhs_; }
  const DiscreteFunctionType &rhs() const { return rhs_; }

  EllipticOperatorType &implicitOperator() { return implicitOperator_; }
  LinearOperatorType &linearOperator() { return linearOperator_; }

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

protected:
  const ModelType &implicitModel_; // the mathematical model

  GridPartType &gridPart_; // grid part(view), e.g. here the leaf grid the discrete space is build with

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_; // the unknown
  DiscreteFunctionType rhs_; // the right hand side

  EllipticOperatorType implicitOperator_; // the implicit operator

  LinearOperatorType linearOperator_; // the linear operator (i.e. jacobian of the implicit)
};

template< class BulkModel, class SurfaceModel, class CoupledGrid >
class CoupledScheme
{
public:
  typedef BulkModel BulkModelType;
  typedef SurfaceModel SurfaceModelType;
  typedef CoupledGrid CoupledGridType;

  typedef typename BulkModelType :: GridPartType BulkGridPartType;
  typedef typename SurfaceModelType :: GridPartType SurfaceGridPartType;

  typedef FemSchemeHolder< BulkModelType, 0 > BulkFemSchemeHolderType;
  typedef FemSchemeHolder< SurfaceModelType, 1 > SurfaceFemSchemeHolderType;

  typedef typename BulkFemSchemeHolderType :: DiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename SurfaceFemSchemeHolderType :: DiscreteFunctionType SurfaceDiscreteFunctionType;

  CoupledScheme( BulkGridPartType &bulkGridPart, SurfaceGridPartType &surfaceGridPart,
		 const BulkModelType &bulkModel, const SurfaceModelType &surfaceModel,
		 const CoupledGridType &coupledGrid )
    : bulkScheme_( bulkGridPart, bulkModel ),
      surfaceScheme_( surfaceGridPart, surfaceModel ),
      coupledGrid_( coupledGrid ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) )
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
    bulk().rhs().clear();
    surface().rhs().clear();

    assembleRHS( bulk().model(), surface().model(),
		 bulkSolution(), surfaceSolution(),
		 coupledGrid_,
		 bulk().rhs(), surface().rhs() );
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

	if( Dune::Fem::MPIManager::rank() == 0 )
	  std::cout << "it: " << iterations_ << " update: " << update
		    << " bulk solver it: " << bulkInvOp.iterations()
		    << " surface solver it: " << surfaceInvOp.iterations()
		    << std::endl;

	prepare();
	++iterations_;
      }
    while( update > eps );
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
  unsigned int iterations_;
};

#endif // #ifndef COUPLED_FEMSCHEME_HH

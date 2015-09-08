#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

// include norms
#include <dune/fem/misc/l2norm.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include header of adaptive scheme
#include "poisson.hh"
#include "coupledscheme.hh"
#include "model.hh"

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class CoupledGridType >
double algorithm ( CoupledGridType &coupledGrid, int step )
{
  // create host grid part consisting of leaf level elements
  typedef typename CoupledGridType :: HBulkGridType HBulkGridType;
  typedef Dune::Fem::AdaptiveLeafGridPart< HBulkGridType, Dune::InteriorBorder_Partition > BulkGridPartType;
  BulkGridPartType bulkGridPart( coupledGrid.bulkGrid() );
  typedef typename CoupledGridType :: HSurfaceGridType HSurfaceGridType;
  typedef Dune::Fem::AdaptiveLeafGridPart< HSurfaceGridType, Dune::InteriorBorder_Partition > SurfaceGridPartType;
  SurfaceGridPartType surfaceGridPart( coupledGrid.surfaceGrid() );

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, HBulkGridType::dimensionworld, 1 > FunctionSpaceType;

  // choose problem
  typedef BulkProblem< FunctionSpaceType > BulkProblemType;
  BulkProblemType bulkProblem;
  typedef SurfaceProblem< FunctionSpaceType > SurfaceProblemType;
  SurfaceProblemType surfaceProblem;

  // type of the mathematical model used
  typedef DiffusionModel< FunctionSpaceType, BulkGridPartType > BulkModelType;
  BulkModelType bulkImplicitModel( bulkProblem, bulkGridPart );
  typedef DiffusionModel< FunctionSpaceType, SurfaceGridPartType > SurfaceModelType;
  SurfaceModelType surfaceImplicitModel( surfaceProblem, surfaceGridPart );

  // create adaptive scheme
  typedef CoupledScheme< BulkModelType, SurfaceModelType, CoupledGridType > SchemeType;
  SchemeType scheme( bulkGridPart, surfaceGridPart, bulkImplicitModel, surfaceImplicitModel, coupledGrid );

  typedef Dune::Fem::GridFunctionAdapter< BulkProblemType, BulkGridPartType > BulkGridExactSolutionType;
  BulkGridExactSolutionType bulkGridExactSolution("bulk exact solution", bulkProblem, bulkGridPart, 5 );
  typedef Dune::Fem::GridFunctionAdapter< SurfaceProblemType, SurfaceGridPartType > SurfaceGridExactSolutionType;
  SurfaceGridExactSolutionType surfaceGridExactSolution("surface exact solution", surfaceProblem, surfaceGridPart, 5 );

  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename SchemeType::BulkDiscreteFunctionType *, BulkGridExactSolutionType * > BulkIOTupleType;
  typedef Dune::Fem::DataOutput< HBulkGridType, BulkIOTupleType > BulkDataOutputType;
  BulkIOTupleType bulkIoTuple( &(scheme.bulkSolution()), &bulkGridExactSolution) ; // tuple with pointers
  BulkDataOutputType bulkDataOutput( coupledGrid.bulkGrid(), bulkIoTuple, DataOutputParameters( step, "bulk" ) );

  typedef Dune::tuple< const typename SchemeType::SurfaceDiscreteFunctionType *, SurfaceGridExactSolutionType * > SurfaceIOTupleType;
  typedef Dune::Fem::DataOutput< HSurfaceGridType, SurfaceIOTupleType > SurfaceDataOutputType;
  SurfaceIOTupleType surfaceIoTuple( &(scheme.surfaceSolution()), &surfaceGridExactSolution) ; // tuple with pointers
  SurfaceDataOutputType surfaceDataOutput( coupledGrid.surfaceGrid(), surfaceIoTuple, DataOutputParameters( step, "surface" ) );

  // setup the right hand side
  scheme.prepare();
  // solve once
  scheme.solve( true );

  bulkDataOutput.write();
  surfaceDataOutput.write();

  // output important info
  std::cout << "h: " << EvolvingDomain :: GridWidth :: gridWidth( bulkGridPart ); << std::endl;
  std::cout << "total dofs: "
	    << scheme.dofs() << std::endl;

  // compute error
  typename Dune :: Fem :: L2Norm< BulkGridPartType > bulkNorm( bulkGridPart );
  const double bulkError = bulkNorm.distance( bulkGridExactSolution, scheme.bulkSolution() );

  typename Dune :: Fem :: L2Norm< SurfaceGridPartType > surfaceNorm( surfaceGridPart );
  const double surfaceError = surfaceNorm.distance( surfaceGridExactSolution, scheme.surfaceSolution() );

  std::cout << "bulkError: " << bulkError << std::endl;
  std::cout << "surfaceError: " << surfaceError << std::endl;

  return std::sqrt( bulkError * bulkError + surfaceError * surfaceError );
}

template< class HBulkGrid, class HSurfaceGrid >
struct CoupledGrid
{
  typedef HBulkGrid HBulkGridType;
  typedef HSurfaceGrid HSurfaceGridType;

  typedef typename HBulkGridType :: LeafGridView BulkLeafGridViewType;

  // type of iterators
  typedef typename BulkLeafGridViewType :: template Codim < 0 > :: Iterator ElementIteratorType;
  typedef typename ElementIteratorType :: Entity BulkEntityType;
  typedef typename BulkEntityType :: EntitySeed BulkEntitySeedType;
  typedef typename BulkEntityType :: EntityPointer BulkEntityPointerType;
  typedef typename BulkLeafGridViewType :: IntersectionIterator IntersectionIteratorType;

  typedef typename Dune :: GridFactory< HSurfaceGridType > GridFactoryType;
  typedef typename std::vector< BulkEntitySeedType > BulkSeedVectorType;

  CoupledGrid( HBulkGrid &bulkGrid )
    : bulkGrid_( bulkGrid )
  {
    assert( HBulkGrid::dimension == 3 );

    // create vertex list to avoid duplicated
    std::vector< Dune :: FieldVector< double, HBulkGrid::dimension > > vertexList;

    // find view of grid
    BulkLeafGridViewType leafView = bulkGrid.leafGridView();

    // add surface grid to factory
    unsigned int count = 0;
    const ElementIteratorType end = leafView.template end< 0 >();
    for( ElementIteratorType it = leafView.template begin< 0 >(); it !=  end;
	 ++ it )   // iterate over elements
      {
	const IntersectionIteratorType iitend = leafView.iend( *it );
	for( IntersectionIteratorType iit = leafView.ibegin( *it ); iit != iitend;
	     ++ iit ) // iterate over intersections
	  {
	    if( iit->boundary() ) // if interseciton is on boundary
	      {
		// add vertices
		std::vector< unsigned int > idx;
		for( int co = 0; co < iit->geometry().corners(); ++co )
		  {
		    // get vertex
		    const Dune :: FieldVector< double, HBulkGrid::dimension >
		      &vertex = iit->geometry().corner( co );

		    // check if already added and find index
		    bool added = false;
		    unsigned int j = 0;
		    for( j = 0; j < vertexList.size(); ++ j )
		      {
			if( ( vertex - vertexList[j] ).two_norm() < 1.0e-8 )
			  {
			    added = true;
			    break;
			  }
		      }

		    if( added )
		      {
			idx.push_back( j );
		      }
		    else
		      {
			factory_.insertVertex( vertex );
			vertexList.push_back( vertex );
			idx.push_back( count );
			++count;
		      }
		  }

		// check if orientation of simplex is correct
		Dune :: FieldVector< double, HBulkGrid::dimension > normal;
		Dune :: FieldVector< double, HBulkGrid::dimension > e1 = vertexList[ idx[0] ] - vertexList[ idx[1] ];
		Dune :: FieldVector< double, HBulkGrid::dimension > e2 = vertexList[ idx[0] ] - vertexList[ idx[2] ];
		normal[0] = e1[1]*e2[2]-e1[2]*e2[1];
		normal[1] = e1[2]*e2[0]-e1[0]*e2[2];
		normal[2] = e1[0]*e2[1]-e1[1]*e2[0];

		if( iit->centerUnitOuterNormal()*normal < 0)
		  std::swap( idx[1], idx[2] );

		// add element to factory
		factory_.insertElement( iit->type(), idx );

		// add bulk element seed to vector
		const BulkEntitySeedType &bulkSeed = it->seed();
		bulkSeedVector_.push_back( bulkSeed );
	      }
	  }
      }

    surfaceGridPtr_ = factory_.createGrid();
  }

  template< class SurfaceEntityType >
  BulkEntitySeedType surfaceBulkMap( const SurfaceEntityType &entity ) const
  {
    const unsigned int seedIndex = factory_.insertionIndex( entity );
    return bulkSeedVector_[ seedIndex ];
  }

  HBulkGridType &bulkGrid() { return bulkGrid_; }
  HSurfaceGridType &surfaceGrid() { return *surfaceGridPtr_; }

private:
  HBulkGridType &bulkGrid_;
  Dune::GridPtr< HSurfaceGridType > surfaceGridPtr_;

  GridFactoryType factory_;
  BulkSeedVectorType bulkSeedVector_;
};

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file
  Dune::Fem::Parameter::append( "../data/parameter" );

  // type of hierarchical grid
  //typedef Dune :: ALUGrid< 3, 3, Dune::simplex, Dune::conforming > HBulkGridType;
  //typedef Dune :: ALUGrid< 2, 3, Dune::simplex, Dune::conforming > HSurfaceGridType;
  typedef Dune :: AlbertaGrid< 3, 3 > HBulkGridType;
  typedef Dune :: AlbertaGrid< 2, 3 > HSurfaceGridType;

  // create grid from DGF file
  const std::string bulkGridkey = Dune::Fem::IOInterface::defaultGridKey( HBulkGridType::dimension );
  const std::string bulkGridFile = Dune::Fem::Parameter::getValue< std::string >( bulkGridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro bulk grid: " << bulkGridFile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HBulkGridType > bulkGridPtr( bulkGridFile );
  HBulkGridType& bulkGrid = *bulkGridPtr ;

  // do initial load balance
  bulkGrid.loadBalance();

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "heat.level" );

  // refine mesh
  const int refineStepsForHalf = Dune::DGFGridInfo< HBulkGridType >::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( bulkGrid, level * refineStepsForHalf );

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "heat.repeats", 0 );

  // calculate first step
  typedef CoupledGrid< HBulkGridType, HSurfaceGridType > CoupledGridType;
  CoupledGridType coupledGrid( bulkGrid );

  double oldError = algorithm( coupledGrid, (repeats > 0) ? 0 : -1 );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( bulkGrid, refineStepsForHalf );
    CoupledGridType coupledGrid( bulkGrid );

    const double newError = algorithm( coupledGrid, step );
    const double eoc = log( oldError / newError ) / M_LN2;
    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      std::cout << "Error: " << newError << std::endl;
      std::cout << "EOC( " << step << " ) = " << eoc << std::endl;
    }
    oldError = newError;
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}

#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

#define DEFORMATION 1
// include geometrty grid part
#include <dune/fem/gridpart/geogridpart.hh>
// include description of surface deformation
#include "deformation.hh"

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include eoc output
#include <dune/fem/misc/femeoc.hh>
#include <dune/evolving-domains/gridwidth.hh>

// include header of adaptive scheme
#include "coupledgrid.hh"
#include "coupledheat.hh"
#include "heatmodel.hh"
#include "coupledheatscheme.hh"

// Assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class CoupledGridType >
void algorithm ( CoupledGridType &coupledGrid, int step, const int eocId )
{
  typedef typename CoupledGridType :: HBulkGridType HBulkGridType;
  typedef typename CoupledGridType :: HSurfaceGridType HSurfaceGridType;

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, HBulkGridType::dimensionworld, 1 > FunctionSpaceType;
  // create time provider
  Dune::Fem::GridTimeProvider< HBulkGridType > timeProvider( coupledGrid.bulkGrid() );

#if DEFORMATION
  // create host grid part consisting of leaf level elements
  typedef Dune::Fem::AdaptiveLeafGridPart< HBulkGridType, Dune::InteriorBorder_Partition > BulkHostGridPartType;
  BulkHostGridPartType bulkHostGridPart( coupledGrid.bulkGrid() );
  typedef Dune::Fem::AdaptiveLeafGridPart< HSurfaceGridType, Dune::InteriorBorder_Partition > SurfaceHostGridPartType;
  SurfaceHostGridPartType surfaceHostGridPart( coupledGrid.surfaceGrid() );

  // construct deformaiton
  typedef DeformationCoordFunction< HBulkGridType :: dimensionworld > DeformationType;
  DeformationType deformation;

  // bulk geometry grid part type
  typedef Dune :: Fem :: GridFunctionAdapter< DeformationType, BulkHostGridPartType > BulkDiscreteDeformationType;
  typedef Dune :: Fem :: GeoGridPart< BulkDiscreteDeformationType > BulkGridPartType;
  BulkDiscreteDeformationType bulkDiscreteDeformation( "bulk deformation", deformation, bulkHostGridPart, 1 );
  BulkGridPartType bulkGridPart( bulkDiscreteDeformation );

  // surface geometry grid part type
  typedef Dune :: Fem :: GridFunctionAdapter< DeformationType, SurfaceHostGridPartType > SurfaceDiscreteDeformationType;
  typedef Dune :: Fem :: GeoGridPart< SurfaceDiscreteDeformationType > SurfaceGridPartType;
  SurfaceDiscreteDeformationType surfaceDiscreteDeformation( "surface deformation", deformation, surfaceHostGridPart, 1 );
  SurfaceGridPartType surfaceGridPart( surfaceDiscreteDeformation );

  typedef CoupledGeoGridPart< CoupledGridType, BulkGridPartType, SurfaceGridPartType > CoupledGeoGridPartType;
  CoupledGeoGridPartType coupledGeoGridPart( coupledGrid, bulkGridPart, surfaceGridPart );
#else
  // create host grid part consisting of leaf level elements
  typedef Dune::Fem::AdaptiveLeafGridPart< HBulkGridType, Dune::InteriorBorder_Partition > BulkGridPartType;
  BulkGridPartType bulkGridPart( coupledGrid.bulkGrid() );
  typedef Dune::Fem::AdaptiveLeafGridPart< HSurfaceGridType, Dune::InteriorBorder_Partition > SurfaceGridPartType;
  SurfaceGridPartType surfaceGridPart( coupledGrid.surfaceGrid() );
#endif

  // type of the mathematical model used
  using BulkHeatModelType = HeatModel< FunctionSpaceType, BulkGridPartType >;
  using SurfaceHeatModelType = HeatModel< FunctionSpaceType, SurfaceGridPartType >;

  // choose problem
  using BulkProblemType = typename BulkHeatModelType :: ProblemType;
  BulkProblemType* bulkProblemPtr = 0;
  using SurfaceProblemType = typename SurfaceHeatModelType :: ProblemType;
  SurfaceProblemType* surfaceProblemPtr = 0;

  const std::string problemNames [] = { "coupled_heat", "coupled_parabolic" };
  const int problemNumber = Dune :: Fem :: Parameter :: getEnum( "coupled.problem", problemNames );
  switch( problemNumber )
    {
    case 0:
      bulkProblemPtr = new BulkHeatProblem< FunctionSpaceType >( timeProvider );
      surfaceProblemPtr = new SurfaceHeatProblem< FunctionSpaceType >( timeProvider );
      break;
    case 1:
      bulkProblemPtr = new BulkParabolicProblem< FunctionSpaceType >( timeProvider );
      surfaceProblemPtr = new SurfaceParabolicProblem< FunctionSpaceType >( timeProvider );
      break;
    default:
      std::cerr << "unrecognised problem name" << std::endl;
    }

  // recover problems
  assert( bulkProblemPtr );
  assert( surfaceProblemPtr );
  BulkProblemType& bulkProblem = *bulkProblemPtr;
  SurfaceProblemType& surfaceProblem = *surfaceProblemPtr;
  typedef ExchangeHeatProblem< FunctionSpaceType > ExchangeProblemType;
  ExchangeProblemType exchangeProblem( timeProvider );

  BulkHeatModelType bulkImplicitModel( bulkProblem, bulkGridPart, true );
  BulkHeatModelType bulkExplicitModel( bulkProblem, bulkGridPart, false );
  SurfaceHeatModelType surfaceImplicitModel( surfaceProblem, surfaceGridPart, true );
  SurfaceHeatModelType surfaceExplicitModel( surfaceProblem, surfaceGridPart, false );
  SurfaceHeatModelType exchangeImplicitModel( exchangeProblem, surfaceGridPart, true );

  // create adaptive scheme
  typedef CoupledHeatScheme< BulkHeatModelType, BulkHeatModelType,
			     SurfaceHeatModelType, SurfaceHeatModelType,
			     SurfaceHeatModelType,
			     CoupledGeoGridPartType > SchemeType;

  typename SchemeType :: BulkFemSchemeHolderType bulkScheme( bulkGridPart, bulkImplicitModel, bulkExplicitModel );
  typename SchemeType :: SurfaceFemSchemeHolderType surfaceScheme( surfaceGridPart, surfaceImplicitModel, surfaceExplicitModel );

  SchemeType scheme( bulkScheme, surfaceScheme,
		     exchangeImplicitModel, coupledGeoGridPart, step );

  typedef Dune::Fem::GridFunctionAdapter< BulkProblemType, BulkGridPartType > BulkGridExactSolutionType;
  BulkGridExactSolutionType bulkGridExactSolution("bulk exact solution", bulkProblem, bulkGridPart, 5 );
  typedef Dune::Fem::GridFunctionAdapter< SurfaceProblemType, SurfaceGridPartType > SurfaceGridExactSolutionType;
  SurfaceGridExactSolutionType surfaceGridExactSolution("surface exact solution", surfaceProblem, surfaceGridPart, 5 );

  //! input/output tuple and setup data writter
  typedef Dune::tuple< const typename SchemeType::BulkDiscreteFunctionType *, BulkGridExactSolutionType * > BulkIOTupleType;
  typedef Dune::Fem::DataOutput< HBulkGridType, BulkIOTupleType > BulkDataOutputType;
  BulkIOTupleType bulkIoTuple( &(scheme.bulkSolution()), &bulkGridExactSolution) ; // tuple with pointers
  BulkDataOutputType bulkDataOutput( coupledGrid.bulkGrid(), bulkIoTuple, DataOutputParameters( step, "bulk" ) );

  typedef Dune::tuple< const typename SchemeType::SurfaceDiscreteFunctionType *, SurfaceGridExactSolutionType * > SurfaceIOTupleType;
  typedef Dune::Fem::DataOutput< HSurfaceGridType, SurfaceIOTupleType > SurfaceDataOutputType;
  SurfaceIOTupleType surfaceIoTuple( &(scheme.surfaceSolution()), &surfaceGridExactSolution) ; // tuple with pointers
  SurfaceDataOutputType surfaceDataOutput( coupledGrid.surfaceGrid(), surfaceIoTuple, DataOutputParameters( step, "surface" ) );

  const double endTime  = Dune::Fem::Parameter::getValue< double >( "heat.endtime", 2.0 );
  const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("heat.reducetimestepfactor", 0.25 );
  double timeStep = Dune::Fem::Parameter::getValue< double >( "heat.timestep", 0.125 );

  timeStep *= pow(dtreducefactor,step);

  //! [time loop]
  // initialize with fixed time step
  timeProvider.init( timeStep ) ;

  // initialize scheme and output initial data
  scheme.initialize();
  // write initial solve
  bulkDataOutput.write( timeProvider );
  surfaceDataOutput.write( timeProvider );
  scheme.closeTimestep( bulkGridExactSolution, surfaceGridExactSolution, timeProvider, true );

  // time loop, increment with fixed time step
  for( ; timeProvider.time() < endTime; timeProvider.next( timeStep ) )
  //! [time loop]
  {
    // assemble explicit pare
    scheme.prepare();

    //! [Set the new time to move to new surface]
    deformation.setTime( timeProvider.time() + timeProvider.deltaT() );

    // solve once - but now we need to reassmble
    scheme.solve(true);

    //! [Set the new time to move to new surface]
    bulkDataOutput.write( timeProvider );
    surfaceDataOutput.write( timeProvider );

    // finalise (compute errors)
    scheme.closeTimestep( bulkGridExactSolution, surfaceGridExactSolution, timeProvider );
  }

  // output final solution
  bulkDataOutput.write( timeProvider );
  surfaceDataOutput.write( timeProvider );


  // get errors
  std::vector< double > store;
  store.push_back( scheme.linftyl2BulkError() );
  store.push_back( scheme.l2h1BulkError() );
  store.push_back( scheme.linftyl2SurfaceError() );
  store.push_back( scheme.l2h1SurfaceError() );
  Dune :: Fem :: FemEoc :: setErrors( eocId, store );

  // write to file / output
  const double h = EvolvingDomain :: GridWidth :: gridWidth( bulkGridPart );
  const int dofs = scheme.nDofs();
  Dune::Fem::FemEoc::write( h, dofs, 0.0, 0.0, std::cout );
}

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

  // initialzie eoc file
  std::string eocOutPath = Dune::Fem::Parameter::getValue<std::string>("fem.eocOutPath", std::string("."));
  Dune::Fem::FemEoc::initialize( eocOutPath, "eoc", "surface only" );

  // add entries to eoc calculation
  std::vector<std::string> femEocHeaders;
  femEocHeaders.push_back("$L^\\infty( L^2(\\Omega) )$ error");
  femEocHeaders.push_back("$L^2( H^1(\\Omega) )$ error");
  femEocHeaders.push_back("$L^\\infty( L^2(\\Gamma) )$ error");
  femEocHeaders.push_back("$L^2( H^1(\\Gamma) )$ error");
  const int eocId = Dune::Fem::FemEoc::addEntry( femEocHeaders );

  // type of hierarchical grid
  //typedef Dune :: ALUGrid< 3, 3, Dune::simplex, Dune::conforming > HBulkGridType;
  //typedef Dune :: ALUGrid< 2, 3, Dune::simplex, Dune::conforming > HSurfaceGridType;
  typedef Dune :: AlbertaGrid< GRIDDIM, WORLDDIM > HBulkGridType;
  typedef Dune :: AlbertaGrid< GRIDDIM-1, WORLDDIM > HSurfaceGridType;

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

  algorithm( coupledGrid, (repeats > 0) ? 0 : -1, eocId );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( bulkGrid, refineStepsForHalf );
    CoupledGridType coupledGrid( bulkGrid );

    algorithm( coupledGrid, step, eocId );
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}

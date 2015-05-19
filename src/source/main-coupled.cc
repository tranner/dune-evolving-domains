#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/filteredgridpart.hh>
#include "../../common/source/radialfilter.hh"
#if FEM_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#else
#include <dune/fem/gridpart/leafgridpart.hh>
#endif

#include <dune/fem_fem_coupling/surfacegridclass.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include header of adaptive scheme
#include "../../common/source/poisson.hh"
#include "../../common/source/femscheme.hh"

#include <dune/fem_fem_coupling/fem_fem_coupling.hh>

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
double algorithm ( HGridType &grid, int step )
{
  // create host grid part consisting of leaf level elements
  typedef typename Dune::FemFemCoupling::Glue<HGridType>::HostGridPartType HostGridPartType;
  HostGridPartType hostGridPart( grid );

  // create filters
  typedef typename Dune::FemFemCoupling::Glue<HGridType>::FilterType FilterType;
  typename FilterType::GlobalCoordinateType center( 0.5 );
  // dummey filter just captures whole of mesh given in
  typename FilterType::ctype radias( 10.0 );
  FilterType filter( hostGridPart, center, radias, true );
  // true filter to generate mesh to be used in glue object generation for overlapping case
  // typename FilterType::ctype radius( -SliceStore::sliceHome() );
  // FilterType Filter( hostGridPart, center, radius, true );

  // we want to solve the problem on the leaf elements of the grid
  typedef typename Dune::FemFemCoupling::Glue< HGridType >::GridPartType GridPartType;
  GridPartType gridPart( hostGridPart, filter );
  // GridPartType GridPart( hostGridPart, Filter );

  // use a scalar function space
  typedef typename Dune::FemFemCoupling::Glue<HGridType>::FunctionSpaceType FunctionSpaceType;

  // type of the mathematical model used
  typedef DiffusionModel< FunctionSpaceType, GridPartType, FunctionSpaceType > ModelType;

  typedef typename ModelType::ProblemType ProblemType ;
  ProblemType* problemPtr = 0 ;
  const std::string problemNames [] = { "cos", "sphere", "sin", "bulkprob", "surfprob" };
  const int problemNumber = Dune::Fem::Parameter::getEnum("volume.problem", problemNames, 0 );
  switch ( problemNumber )
  {
    case 0:
      problemPtr = new CosinusProduct< FunctionSpaceType > ();
      break ;
    case 1:
      problemPtr = new SphereProblem< FunctionSpaceType > ();
      break ;
    case 2:
      problemPtr = new SinusProduct< FunctionSpaceType > ();
      break ;
    case 3:
      problemPtr = new BulkProblem< FunctionSpaceType > ();
      break ;
    case 4:
      problemPtr = new SurfaceProblem< FunctionSpaceType > ();
      break ;
    default:
      problemPtr = new CosinusProduct< FunctionSpaceType > ();
  }
  assert( problemPtr );
  ProblemType& problem = *problemPtr ;

  // implicit model for left hand side
  ModelType implicitModel( problem, gridPart );
  // ModelType ImplicitModel( problem, GridPart );

  // create adaptive scheme
  typedef FemScheme< ModelType, istl > SchemeType;
  SchemeType scheme( gridPart, implicitModel,"bulk" );

  // create dune grid file with surface mesh
  typedef typename Dune::FemFemCoupling::Glue< HGridType >::GrydType GrydType;
  // create dune grid file with surface mesh
  GrydType* grydPtr = scheme.extractSurface();
  GrydType& gryd = *grydPtr ;

  // do initial load balance
  gryd.loadBalance();
  
  typedef typename Dune::FemFemCoupling::Glue< HGridType >::GrydPartType GrydPartType;
  GrydPartType grydPart(gryd);
  
  // use a scalar function space
  typedef typename Dune::FemFemCoupling::Glue<HGridType>::FunctionSpaceOutsideType ExtractedSurfaceFunctionSpaceType;

  // type of the mathematical model used
  typedef DiffusionModel< ExtractedSurfaceFunctionSpaceType, GrydPartType, ExtractedSurfaceFunctionSpaceType > ExtractedSurfaceModelType;

  typedef typename ExtractedSurfaceModelType::ProblemType ExtractedSurfaceProblemType ;
  ExtractedSurfaceProblemType* problimPtr = 0 ;

  const int problimNumber = Dune::Fem::Parameter::getEnum("surface.problem", problemNames, 0 );

  switch ( problimNumber )
  {
    case 0:
      problimPtr = new CosinusProduct< ExtractedSurfaceFunctionSpaceType > ();
      break ;
    case 1:
      problimPtr = new SphereProblem< ExtractedSurfaceFunctionSpaceType > ();
      break ;
    case 2:
      problimPtr = new SinusProduct< ExtractedSurfaceFunctionSpaceType > ();
      break ;
    case 3:
      problimPtr = new BulkProblem< ExtractedSurfaceFunctionSpaceType > ();
      break ;
    case 4:
      problimPtr = new SurfaceProblem< ExtractedSurfaceFunctionSpaceType > ();
      break ;
    default:
      problimPtr = new CosinusProduct< ExtractedSurfaceFunctionSpaceType > ();
  }
  assert( problemPtr );
  ExtractedSurfaceProblemType& problim = *problimPtr ;

  ExtractedSurfaceModelType implycitModel( problim, grydPart );
  
  typedef FemScheme< ExtractedSurfaceModelType, armadillo > ExtractedSurfaceSchemeType;
  ExtractedSurfaceSchemeType schyme( grydPart, implycitModel,"surface" );

  // SchemeType scheme( gridPart, implicitModel );
  // int ghg = 0;
  // std::cin >> ghg;

  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GrydPartType > GrydExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
  GrydExactSolutionType grydExactSolution("exact solution", problim, grydPart, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution) ; // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  typedef Dune::tuple< const typename ExtractedSurfaceSchemeType::DiscreteFunctionType *, GrydExactSolutionType *  > SurfaceIOTupleType;
  typedef Dune::Fem::DataOutput< GrydType, SurfaceIOTupleType > SurfaceDataOutputType;
  SurfaceIOTupleType surfaceioTuple( &(schyme.solution()), &grydExactSolution ) ; // tuple with pointers
  SurfaceDataOutputType surfacedataOutput( gryd, surfaceioTuple, DataOutputParameters( step, "surface" ) );

  // calculate error 
  double error = 0;
  double errer = 0;

#if 0
  // iterate a certain number of times
  for( int n = 0; n <= 31; ++n )
  {
    // set-up the coupling
    scheme.setup(schyme.solution());
    schyme.setup(scheme.solution());

    // setup the right hand side
    scheme.prepare();
    schyme.prepare();

    // apply the coupling
    scheme.couple(schyme.solution());
    schyme.couple(scheme.solution());

    // solve once 
    scheme.solve( n == 0 );
    schyme.solve( n == 0 );

    // calculate standard error 
    // select norm for error computation
    typedef Dune::Fem::L2Norm< GridPartType > NormType;
    typedef Dune::Fem::L2Norm< GrydPartType > NurmType;
    NormType norm( gridPart );
    NurmType nurm( grydPart );
    error = norm.distance( gridExactSolution, scheme.solution() );
    errer = nurm.distance( grydExactSolution, schyme.solution() );
    std::cout << "                                                   Error at iteration " << n << " = " << error << " + " << errer << " = " << error + errer << std::endl;

    // check for convergence
    bool done = scheme.stop(); bool dune = schyme.stop();
    if ( ( n > 3 ) && done && dune )
    {
      std::cout << "Converged at iteration " << n << " to an error of " << error+errer << std::endl;
      break;
    }

    // only write output (and continue iterating) if solution still converging
    else
    {
      // write initial solve 
      dataOutput.write();
      surfacedataOutput.write();

      // reset the problem
      scheme.reset();
      schyme.reset();
    }
  }
#endif
  // write initial solve 
  dataOutput.write();
  surfacedataOutput.write();

  return error + errer;
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

  // type of hierarchical grid
  typedef Dune::GridSelector::GridType  HGridType ;

  // set slice value to be used for the gluing
  Dune::SliceStore::sliceHome(-Dune::Fem::Parameter::getValue< double >( "coupling.inner" ));
  Dune::SliceStore::x(1); Dune::SliceStore::y(1); 
  if( HGridType::dimension < 3 )
  {
    Dune::SliceStore::z(0);
  }
  else
  {
    Dune::SliceStore::z(1);
  }
  Dune::SliceStore::sliceHome(-Dune::Fem::Parameter::getValue< double >( "coupling.outer" ),true);

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "poisson.level" );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "poisson.repeats", 0 );

  // calculate first step
  double oldError = algorithm( grid, (repeats > 0) ? 0 : -1 );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    const double newError = algorithm( grid, step );
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

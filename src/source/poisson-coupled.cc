#include <config.h>

// iostream includes
#include <iostream>

// dune grid includes
#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#define DEFORMATION 1
// include geometrty grid part
#if USE_OLD_GRIDPART
#include <dune/fem/gridpart/geogridpart.hh>
#else
#include <dune/evolving-domains/gridpart/geogridpart.hh>
#endif
// include description of surface deformation
#include <dune/evolving-domains/deformation.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include eoc output
#include <dune/fem/misc/femeoc.hh>
#include <dune/evolving-domains/gridwidth.hh>

// include header of adaptive scheme
#include "coupledgrid.hh"
#include "coupledpoisson.hh"

#include "model.hh"
#include "coupledscheme.hh"

#include <dune/fem/space/common/functionspace.hh>
template< int dimWorld >
struct BoundaryProjection
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
    y /= y.two_norm();
  }

  bool onBoundary( const DomainType& x )
  {
    return std::abs( x.two_norm() - 1.0 ) < 1.0e-12;
  }
};

template< int dimWorld >
struct DeformationCoordFunction
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  explicit DeformationCoordFunction ()
  {}

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
  }
};

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


  // construct deformation
  typedef BoundaryProjection< HBulkGridType::dimensionworld > BoundaryProjectionType;
  BoundaryProjectionType boundaryProjection;
  typedef DeformationCoordFunction< HBulkGridType::dimensionworld > DeformationType;
  DeformationType deformation;

  typedef DiscreteDeformationCoordHolder< DeformationType, BoundaryProjectionType,
					  BulkHostGridPartType, 0, POLORDER > BulkDiscreteDeformationCoordHolderType;
  typedef typename BulkDiscreteDeformationCoordHolderType :: DiscreteFunctionType BulkCoordinateFunctionType;
  BulkDiscreteDeformationCoordHolderType bulkDiscreteDeformation( deformation, boundaryProjection, bulkHostGridPart );

  // surface geometry grid part type
  typedef Dune :: Fem :: GeoGridPart< BulkCoordinateFunctionType > BulkGridPartType;
  BulkGridPartType bulkGridPart( bulkDiscreteDeformation.coordFunction() );

  typedef DiscreteDeformationCoordHolder< DeformationType, BoundaryProjectionType,
					  SurfaceHostGridPartType, 1, POLORDER > SurfaceDiscreteDeformationCoordHolderType;
  typedef typename SurfaceDiscreteDeformationCoordHolderType :: DiscreteFunctionType SurfaceCoordinateFunctionType;
  SurfaceDiscreteDeformationCoordHolderType surfaceDiscreteDeformation( deformation, boundaryProjection, surfaceHostGridPart );

  // surface geometry grid part type
  typedef Dune :: Fem :: GeoGridPart< SurfaceCoordinateFunctionType > SurfaceGridPartType;
  SurfaceGridPartType surfaceGridPart( surfaceDiscreteDeformation.coordFunction() );

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
  using BulkModelType = DiffusionModel< FunctionSpaceType, BulkGridPartType >;
  using SurfaceModelType = DiffusionModel< FunctionSpaceType, SurfaceGridPartType >;

  // choose problem
  using BulkProblemType = typename BulkModelType :: ProblemType;
  BulkProblemType* bulkProblemPtr = 0;
  using SurfaceProblemType = typename SurfaceModelType :: ProblemType;
  SurfaceProblemType* surfaceProblemPtr = 0;

  // const std::string problemNames [] = { "coupled_poisson" };
  // const int problemNumber = Dune :: Fem :: Parameter :: getEnum( "coupledpoisson.problem", problemNames );
  const int problemNumber = 0;
  switch( problemNumber )
    {
    case 0:
      bulkProblemPtr = new BulkProblem< FunctionSpaceType >();
      surfaceProblemPtr = new SurfaceProblem< FunctionSpaceType >();
      break;
    default:
      std::cerr << "unrecognised problem name" << std::endl;
    }

  // recover problems
  assert( bulkProblemPtr );
  assert( surfaceProblemPtr );
  BulkProblemType& bulkProblem = *bulkProblemPtr;
  SurfaceProblemType& surfaceProblem = *surfaceProblemPtr;
  typedef ExchangeProblem< FunctionSpaceType > ExchangeProblemType;
  ExchangeProblemType exchangeProblem;

  BulkModelType bulkImplicitModel( bulkProblem, bulkGridPart );
  SurfaceModelType surfaceImplicitModel( surfaceProblem, surfaceGridPart );
  SurfaceModelType exchangeImplicitModel( exchangeProblem, surfaceGridPart );

  // create adaptive scheme
  typedef CoupledScheme< BulkModelType, SurfaceModelType,
			 SurfaceModelType, CoupledGeoGridPartType > SchemeType;

  typename SchemeType :: BulkFemSchemeHolderType bulkScheme( bulkGridPart, bulkImplicitModel );
  typename SchemeType :: SurfaceFemSchemeHolderType surfaceScheme( surfaceGridPart, surfaceImplicitModel );

  SchemeType scheme( bulkScheme, surfaceScheme,
		     exchangeImplicitModel, coupledGeoGridPart );

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

  // assemble explicit pare
  scheme.prepare();

  // solve once - but now we need to reassmble
  scheme.solve(true);

  // output final solution
  bulkDataOutput.write( timeProvider );
  surfaceDataOutput.write( timeProvider );

  // get errors
  const std::vector< double > store = {
    bulkScheme.l2Error( bulkGridExactSolution ),
    bulkScheme.h1Error( bulkGridExactSolution ),
    surfaceScheme.l2Error( surfaceGridExactSolution ),
    surfaceScheme.h1Error( surfaceGridExactSolution )
  };
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
  femEocHeaders.push_back("$L^2(\\Omega)$ error");
  femEocHeaders.push_back("$H^1(\\Omega)$ error");
  femEocHeaders.push_back("$L^2(\\Gamma)$ error");
  femEocHeaders.push_back("$H^1(\\Gamma)$ error");
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
  const int level = Dune::Fem::Parameter::getValue< int >( "poisson.level" );

  // refine mesh
  const int refineStepsForHalf = Dune::DGFGridInfo< HBulkGridType >::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( bulkGrid, level * refineStepsForHalf );

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "poisson.repeats", 0 );

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

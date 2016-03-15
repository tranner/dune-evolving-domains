#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

// refinement includes
#include <dune/fem/space/common/adaptmanager.hh>

// include geometrty grid part
#include <dune/fem/gridpart/geogridpart.hh>
// include description of surface deformation
#include "deformation.hh"

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include eoc output
#include <dune/fem/misc/femeoc.hh>
#include <dune/evolving-domains/gridwidth.hh>

// quadrature for area calculation
#include <dune/fem/quadrature/cachingquadrature.hh>

#include "coupledgrid.hh"

template< class GridPartType >
void computeArea( const GridPartType& gridPart, double& bulkArea, double& surfaceArea )
{
  // compute area using quadrature rules
  bulkArea = 0.0;
  surfaceArea = 0.0;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef typename GridPartType::GridViewType GridViewType;
  for( const auto& e : elements( static_cast<GridViewType>( gridPart ) ) )
    {
      const auto &geometry = e.geometry();

      QuadratureType quadrature( e, POLORDER+1 );
      const size_t numQuadraturePoints = quadrature.nop();

      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  //! [Compute local contribution of operator]
	  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

	  bulkArea += weight;
	}

      const auto& iitEnd = gridPart.iend( e );
      for( auto iit = gridPart.ibegin( e ); iit != iitEnd; ++iit )
	{
	  const auto& intersection = *iit;
	  if( intersection.boundary() )
	    {
	      // extract intersection geometry
	      const auto& intersectionGeometry = intersection.geometry();

	      // find quadrature rule on intersection
	      FaceQuadratureType quadInside( gridPart, intersection, POLORDER+1,
					     FaceQuadratureType :: INSIDE );
	      const size_t numQuadraturePoints = quadInside.nop();

	      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
		{
		  // find quadrature
		  const typename FaceQuadratureType::LocalCoordinateType &xLocal = quadInside.localPoint( pt );
		  double weight = quadInside.weight( pt ) * intersectionGeometry.integrationElement( xLocal );

		  surfaceArea += weight;
		}
	    }
	}
    }
}

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

  // create host grid part consisting of leaf level elements
  typedef Dune::Fem::AdaptiveLeafGridPart< HBulkGridType, Dune::InteriorBorder_Partition > BulkHostGridPartType;
  BulkHostGridPartType bulkHostGridPart( coupledGrid.bulkGrid() );
  typedef Dune::Fem::AdaptiveLeafGridPart< HSurfaceGridType, Dune::InteriorBorder_Partition > SurfaceHostGridPartType;
  SurfaceHostGridPartType surfaceHostGridPart( coupledGrid.surfaceGrid() );

  // construct deformaiton
  typedef DeformationCoordFunction< HBulkGridType :: dimensionworld > DeformationType;
  DeformationType deformation;

  typedef DiscreteDeformationCoordHolder< DeformationType, BulkHostGridPartType, 1 > DiscreteDeformationCoordHolderType;
  typedef typename DiscreteDeformationCoordHolderType :: DiscreteFunctionType CoordinateFunctionType;
  DiscreteDeformationCoordHolderType discreteDeformation( deformation, bulkHostGridPart );

  // bulk geometry grid part type
  typedef Dune :: Fem :: GeoGridPart< CoordinateFunctionType > BulkGridPartType;
  BulkGridPartType bulkGridPart( discreteDeformation.coordFunction() );

  // compute area
  double bulkArea, surfaceArea;
  computeArea( bulkGridPart, bulkArea, surfaceArea );

  // get errors
  std::vector< double > store = { std::abs( bulkArea - 4.0 * M_PI / 3.0 ),
				  std::abs( surfaceArea - 4.0 * M_PI ) };
  Dune :: Fem :: FemEoc :: setErrors( eocId, store );

  // write to file / output
  const double h = EvolvingDomain :: GridWidth :: gridWidth( bulkGridPart );
  const int dofs = 0;
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
  Dune::Fem::Parameter::append( "parameter" );

  // initialzie eoc file
  std::string eocOutPath = Dune::Fem::Parameter::getValue<std::string>("fem.eocOutPath", std::string("."));
  Dune::Fem::FemEoc::initialize( eocOutPath, "eoc", "surface only" );

  // add entries to eoc calculation
  std::vector<std::string> femEocHeaders { "bulk area error", "boundary area error" };
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
  const int level = Dune::Fem::Parameter::getValue< int >( "areatest.level" );

  // refine mesh
  const int refineStepsForHalf = Dune::DGFGridInfo< HBulkGridType >::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( bulkGrid, level * refineStepsForHalf );

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "areatest.repeats", 0 );

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

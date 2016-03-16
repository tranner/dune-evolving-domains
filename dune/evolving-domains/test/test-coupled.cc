#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// refinement includes
#include <dune/fem/space/common/adaptmanager.hh>

// include geometrty grid part
#if USE_OLD_GRIDPART
#include <dune/fem/gridpart/geogridpart.hh>
#else
#include <dune/evolving-domains/gridpart/geogridpart.hh>
#endif
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

// timer
#include <ctime>

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
struct Deformation
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
  }
};

template< class GridPart >
double computeArea( const GridPart& gridPart )
{
  // compute area using quadrature rules
  double ret = 0;
  typedef Dune::Fem::CachingQuadrature< GridPart, 0 > QuadratureType;
  typedef typename GridPart::GridViewType GridViewType;
  for( const auto& e : elements( static_cast<GridViewType>( gridPart ) ) )
    {
      const auto &geometry = e.geometry();

      QuadratureType quadrature( e, 2*POLORDER+1 );
      const size_t numQuadraturePoints = quadrature.nop();

      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  //! [Compute local contribution of operator]
	  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

	  ret += weight;
	}
    }

  return ret;
}

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

  // construct deformation
  typedef BoundaryProjection< HBulkGridType::dimensionworld > BoundaryProjectionType;
  BoundaryProjectionType boundaryProjection;
  typedef Deformation< HBulkGridType::dimensionworld > DeformationType;
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

  // start timer
  const std::clock_t start = std::clock();

  // compute area
  double bulkArea, boundaryArea;
  computeArea( bulkGridPart, bulkArea, boundaryArea );

  const double surfaceArea = computeArea( surfaceGridPart );

  const double timeEllapsed = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

  // get errors
  std::vector< double > store = { std::abs( bulkArea - 4.0 * M_PI / 3.0 ),
				  std::abs( boundaryArea - 4.0 * M_PI ),
				  std::abs( surfaceArea - 4.0 * M_PI ) };
  std::cout << "bulk area:\t" << std::setprecision (15) << bulkArea << std::endl;
  std::cout << "boundary area:\t" << std::setprecision (15) << boundaryArea << std::endl;
  std::cout << "surface area:\t" << std::setprecision (15) << surfaceArea << std::endl;
  Dune :: Fem :: FemEoc :: setErrors( eocId, store );

  // write to file / output
  const double h = EvolvingDomain :: GridWidth :: gridWidth( bulkGridPart );
  const int dofs = bulkDiscreteDeformation.dofs();
  const int elements = bulkDiscreteDeformation.elements();
  Dune::Fem::FemEoc::write( h, dofs, timeEllapsed, elements, std::cout );
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
  std::vector<std::string> femEocHeaders { "bulk area error", "boundary area error", "surface area error" };
  const int eocId = Dune::Fem::FemEoc::addEntry( femEocHeaders );

  // type of hierarchical grid
  typedef Dune :: GridSelector :: GridType HBulkGridType;
#ifdef ALBERTAGRID
  typedef Dune :: AlbertaGrid< GRIDDIM-1, WORLDDIM > HSurfaceGridType;
#endif
#ifdef ALUGRID_CONFORM
  typedef Dune :: ALUGrid< GRIDDIM-1, WORLDDIM, Dune::simplex, Dune::conforming > HSurfaceGridType;
#endif

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

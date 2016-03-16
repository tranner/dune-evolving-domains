#include <config.h>

// iostream includes
#include <iostream>

// grid includes
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>

#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

template <class HGridType >
void algorithm ( HGridType &grid, const int step )
{
  int n = 0;
  double volume = 0;

  const auto end = grid.template leafend<0>();
  for( auto it = grid.template leafbegin<0>(); it != end; ++it )
    {
      volume += it->geometry().volume();
      n++;
    }

  std::cout << "level: " << step
	    << " elements: " << n
	    << " volume: " << volume
	    << " error: " << std::abs( volume - 4.0 * M_PI / 3.0 )
	    << std::endl;
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // create grid from DGF file
  const std::string gridFile = "ball.dgf";

  {
    // type of hierarchical grid
    typedef Dune :: ALUGrid< 3, 3, Dune::simplex, Dune::conforming > HGridType;
    std::cout << "Dune :: ALUGrid< 3, 3, Dune::simplex, Dune::conforming >" << std::endl;

    // the method rank and size from MPIManager are static
    std::cout << "Loading macro bulk grid: " << gridFile << std::endl;

    // construct macro using the DGF Parser
    Dune::GridPtr< HGridType > gridPtr( gridFile );
    HGridType& grid = *gridPtr ;

    // do initial load balance
    grid.loadBalance();

    const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

    for( int step = 0; step <= 5; ++step )
      {
	// refine globally such that grid with is bisected
	// and all memory is adjusted correctly
	grid.globalRefine( refineStepsForHalf );

	algorithm( grid, step );
      }
  }

  {
    typedef Dune::AlbertaGrid< 3 > HGridType;
    std::cout << "Dune :: AlbertaGrid< 3 >" << std::endl;

    // the method rank and size from MPIManager are static
    std::cout << "Loading macro bulk grid: " << gridFile << std::endl;

    // construct macro using the DGF Parser
    Dune::GridPtr< HGridType > gridPtr( gridFile );
    HGridType& grid = *gridPtr ;

    // do initial load balance
    grid.loadBalance();

    const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

    for( int step = 0; step <= 5; ++step )
      {
	// refine globally such that grid with is bisected
	// and all memory is adjusted correctly
	grid.globalRefine( refineStepsForHalf );

	algorithm( grid, step );
      }
  }


  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}

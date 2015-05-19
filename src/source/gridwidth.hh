#ifndef GRIDWIDTH_HH
#define GRIDWIDTH_HH

#include <dune/common/binaryfunctions.hh>

namespace EvolvingDomain
{
  // GridWidth
  // ---------
  namespace GridWidth
  {
    /**
     * \function gridWidth
     *
     * \brief calculates maximum element diameter
     */
    template < class GridPart >
    double gridWidth( const GridPart &gridPart )
    {
      return calcGridWidth( gridPart, Dune :: Max<double>() );
    }

    /**
     * \function minGridWidth
     *
     * \brief calculates minimum element diameter
     */
    template < class GridPart >
    double minGridWidth( const GridPart &gridPart )
    {
      return calcGridWidth( gridPart, Dune :: Min<double>() );
    }

    /**
     * \function gridWidth
     *
     * \brief calculates maximum element diameter in zero level set
     */
    template < class GridPart, class LevelSetFunction >
    double gridWidth( const GridPart &gridPart, const LevelSetFunction &levelSetFunction )
    {
      return calcGridWidth( gridPart, Dune :: Max<double>(), levelSetFunction );
    }

    /**
     * \function minGridWidth
     *
     * \brief calculates minimum element diameter in zero level set
     */
    template < class GridPart, class LevelSetFunction >
    double minGridWidth( const GridPart &gridPart, const LevelSetFunction &levelSetFunction )
    {
      return calcGridWidth( gridPart, Dune :: Min<double>(), levelSetFunction );
    }

    template < class MinMax >
    struct MinMaxInit;

    template < class T >
    struct MinMaxInit< Dune :: Min< T > >
    {
      static T init()
      {
	return std::numeric_limits< T >::max();
      }
    };

    template < class T >
    struct MinMaxInit< Dune :: Max< T > >
    {
      static T init()
      {
	return std::numeric_limits< T >::min();
      }
    };

    template < class GridPart, class MinMax >
    double calcGridWidth( const GridPart &gridPart, const MinMax& minmax )
    {
      typedef GridPart GridPartType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;
      typedef typename EntityType::Geometry GeometryType;

      typedef typename GeometryType :: GlobalCoordinate CoordType;

      double h = MinMaxInit< MinMax > :: init();

      const IteratorType end = gridPart.template end<0> ();
      for( IteratorType it = gridPart.template begin<0> (); it != end; ++it )
	{
	  const EntityType &entity = *it;
	  const GeometryType &geometry = entity.geometry();

	  const int nCorners = geometry.corners();

	  for( int i = 0; i < nCorners; ++i )
	    {
	      CoordType corneri = geometry.corner( i );

	      for( int j = 0; j < i; ++j )
		{
		  CoordType cornerj = geometry.corner( j );
		  h = minmax( ( corneri - cornerj ).two_norm(), h );
		}
	    }
	}

      double h_ = h;
      gridPart.grid().comm().template allreduce<MinMax> ( &h_, &h, 1 );

      return h;
    }

    template < class GridPart, class LevelSetFunction, class MinMax >
    double calcGridWidth( const GridPart &gridPart, const MinMax &minmax, const LevelSetFunction &levelSetFunction )
    {
      typedef GridPart GridPartType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;
      typedef typename EntityType::Geometry GeometryType;

      typedef typename GeometryType :: GlobalCoordinate CoordType;

      double h = MinMaxInit< MinMax > :: init();

      const IteratorType end = gridPart.template end<0> ();
      for( IteratorType it = gridPart.template begin<0> (); it != end; ++it )
	{
	  const EntityType &entity = *it;

	  if( levelSetFunction.contains( entity ) )
	    {
	      const GeometryType &geometry = entity.geometry();

	      const int nCorners = geometry.corners();

	      for( int i = 0; i < nCorners; ++i )
		{
		  CoordType corneri = geometry.corner( i );

		  for( int j = 0; j < i; ++j )
		    {
		      CoordType cornerj = geometry.corner( j );
		      h = minmax( ( corneri - cornerj ).two_norm(), h );
		    }
		}
	    }
	}

      double h_ = h;
      gridPart.grid().comm().template allreduce<MinMax> ( &h_, &h, 1 );

      return h;
    }

  }
}

#endif

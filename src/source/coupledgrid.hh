#ifndef COUPLEDGRID_HH
#define COUPLEDGRID_HH

// stl include
#include <vector>

// dune common include
#include <dune/common/fvector.hh>

// dune grid include
#include <dune/grid/albertagrid/gridfactory.hh>


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
  typedef typename BulkLeafGridViewType :: IntersectionIterator IntersectionIteratorType;

  typedef typename HSurfaceGridType :: LeafGridView :: template Codim< 0 > :: Iterator :: Entity SurfaceEntityType;

  typedef typename Dune :: GridFactory< HSurfaceGridType > GridFactoryType;
  typedef typename std::vector< BulkEntitySeedType > BulkSeedVectorType;

  CoupledGrid( HBulkGrid &bulkGrid )
    : bulkGrid_( bulkGrid )
  {
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


		if( HBulkGrid::dimension == 3 )
		  {
		    // check if orientation of simplex is correct
		    Dune :: FieldVector< double, HBulkGrid::dimension > normal;
		    Dune :: FieldVector< double, HBulkGrid::dimension > e1 = vertexList[ idx[0] ] - vertexList[ idx[1] ];
		    Dune :: FieldVector< double, HBulkGrid::dimension > e2 = vertexList[ idx[0] ] - vertexList[ idx[2] ];
		    normal[0] = e1[1]*e2[2]-e1[2]*e2[1];
		    normal[1] = e1[2]*e2[0]-e1[0]*e2[2];
		    normal[2] = e1[0]*e2[1]-e1[1]*e2[0];

		    if( iit->centerUnitOuterNormal()*normal < 0)
		      std::swap( idx[1], idx[2] );
		  }

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

  BulkEntitySeedType surfaceBulkMap( const SurfaceEntityType &entity ) const
  {
    const unsigned int seedIndex = factory_.insertionIndex( entity );
    return bulkSeedVector_[ seedIndex ];
  }
  BulkEntitySeedType surfaceBulkMap( const unsigned int idx ) const
  {
    return bulkSeedVector_[ idx ];
  }

  typename BulkEntityType :: Geometry :: LocalCoordinate
  surfaceLocalToBulkLocal( const SurfaceEntityType& entity, const typename SurfaceEntityType :: Geometry :: LocalCoordinate& x ) const
  {
    // find global coordinate
    const typename SurfaceEntityType :: Geometry surfaceGeometry = entity.geometry();
    const typename SurfaceEntityType :: Geometry :: GlobalCoordinate xGlobal = surfaceGeometry.global( x );

    // find bulk entity
    const BulkEntitySeedType seed = surfaceBulkMap( entity );
    const BulkEntityType bulkEntity = bulkGrid().entity( seed );

    // find local coordinate for bulk geometry
    const typename BulkEntityType :: Geometry bulkGeometry = bulkEntity.geometry();
    const typename BulkEntityType :: Geometry :: LocalCoordinate xBulkLocal = bulkGeometry.local( xGlobal );

    return xBulkLocal;
  }

  HBulkGridType &bulkGrid() { return bulkGrid_; }
  HSurfaceGridType &surfaceGrid() { return *surfaceGridPtr_; }

  const HBulkGridType &bulkGrid() const { return bulkGrid_; }
  const HSurfaceGridType &surfaceGrid() const { return *surfaceGridPtr_; }

  unsigned int seedIndex( const SurfaceEntityType &entity ) const
  {
    return factory_.insertionIndex( entity );
  }

  unsigned int maxSeedIndex() const
  {
    return bulkSeedVector_.size();
  }

private:
  HBulkGridType &bulkGrid_;
  Dune::GridPtr< HSurfaceGridType > surfaceGridPtr_;

  GridFactoryType factory_;
  BulkSeedVectorType bulkSeedVector_;
};

template< class CoupledGrid, class BulkGeoGridPart, class SurfaceGeoGridPart >
struct CoupledGeoGridPart
{
  typedef CoupledGrid HostGridType;
  typedef BulkGeoGridPart BulkGeoGridPartType;
  typedef SurfaceGeoGridPart SurfaceGeoGridPartType;

  typedef typename HostGridType :: BulkEntityType BulkHostEntityType;
  typedef typename HostGridType :: BulkEntitySeedType BulkHostEntitySeedType;
  typedef typename HostGridType :: SurfaceEntityType SurfaceHostEntityType;

  typedef typename BulkGeoGridPartType :: template Codim< 0 > :: IteratorType BulkIteratorType;
  typedef typename BulkIteratorType :: Entity BulkEntityType;
  typedef typename BulkEntityType :: EntitySeed BulkEntitySeedType;
  typedef typename SurfaceGeoGridPartType :: template Codim< 0 > :: IteratorType :: Entity SurfaceEntityType;

  typedef std::vector< BulkEntitySeedType > BulkHostGeoMapType;

  CoupledGeoGridPart( const HostGridType &coupledGrid,
		      const BulkGeoGridPartType &bulkGridPart,
		      const SurfaceGeoGridPartType &surfaceGridPart )
    : coupledGrid_( coupledGrid ),
      bulkGridPart_( bulkGridPart ),
      surfaceGridPart_( surfaceGridPart ),
      map_( coupledGrid_.maxSeedIndex() )
  {
    BulkIteratorType end = bulkGridPart.template end< 0 >();
    for( BulkIteratorType it = bulkGridPart.template begin< 0 >(); it != end; ++it )
      {
	const BulkEntityType &entity = *it;
	const BulkHostEntityType &hostEntity = gridEntity( entity );

	for( unsigned int i = 0; i < coupledGrid_.maxSeedIndex(); ++i )
	  {
	    if( map_[ i ].isValid() )
	      continue;

	    const BulkHostEntitySeedType &otherHostSeed = surfaceBulkMap( i );
	    const BulkHostEntityType &otherHostEntity = coupledGrid_.bulkGrid().entity( otherHostSeed );

	    if( hostEntity == otherHostEntity )
	      {
		map_[ i ] = entity.seed();
	      }
	  }
      }
  }
  BulkEntitySeedType surfaceBulkMap( const SurfaceEntityType& entity ) const
  {
    const unsigned int idx = coupledGrid_.seedIndex( gridEntity( entity ) );
    return map_[ idx ];
  }

  typename BulkEntityType :: Geometry :: LocalCoordinate
  surfaceLocalToBulkLocal( const SurfaceEntityType& entity, const typename SurfaceEntityType :: Geometry :: LocalCoordinate& x ) const
  {
    return coupledGrid_.surfaceLocalToBulkLocal( gridEntity( entity ), x );
  }


  const BulkGeoGridPartType &bulkGridPart() const
  {
    return bulkGridPart_;
  }

  const SurfaceGeoGridPartType &surfaceGridPart() const
  {
    return surfaceGridPart_;
  }

protected:
  BulkHostEntitySeedType surfaceBulkMap( const SurfaceHostEntityType &entity ) const
  {
    return coupledGrid_.surfaceBulkMap( entity );
  }
  BulkHostEntitySeedType surfaceBulkMap( const unsigned int idx ) const
  {
    return coupledGrid_.surfaceBulkMap( idx );
  }

private:
  const HostGridType &coupledGrid_;
  const BulkGeoGridPartType &bulkGridPart_;
  const SurfaceGeoGridPartType &surfaceGridPart_;

  BulkHostGeoMapType map_;
};

#endif // #ifndef COUPLEDGRID_HH

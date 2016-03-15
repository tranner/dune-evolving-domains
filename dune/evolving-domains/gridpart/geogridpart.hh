#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_HH

#if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
#error "Experimental grid extensions required for GeoGridPart. Reconfigure with --enable-experimental-grid-extensions to enable GeoGridPart."
#else

#include <dune/grid/common/gridview.hh>

#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/metatwistutility.hh>
#include <dune/evolving-domains/gridpart/geogridpart/capabilities.hh>
#include <dune/evolving-domains/gridpart/geogridpart/datahandle.hh>
#include <dune/evolving-domains/gridpart/geogridpart/entity.hh>
#include <dune/evolving-domains/gridpart/geogridpart/geometry.hh>
#include <dune/evolving-domains/gridpart/geogridpart/intersection.hh>
#include <dune/evolving-domains/gridpart/geogridpart/intersectioniterator.hh>
#include <dune/evolving-domains/gridpart/geogridpart/iterator.hh>
#include <dune/fem/gridpart/idgridpart/indexset.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class CoordFunction >
    class GeoGridPart;

    template< class CoordFunction >
    class GeoGridPartFamily;



    // GeoGridPartFamily
    // -----------------

    // Traits for dune-grid facades ("Gen-Gurke!")
    template< class CoordFunction >
    class GeoGridPartFamily
    {
    public:
      typedef typename CoordFunction::RangeFieldType ctype;

      static const int dimension = CoordFunction::GridPartType::dimension;
      static const int dimensionworld = CoordFunction::FunctionSpaceType::dimRange;

      typedef GeoGridPartFamily< CoordFunction > GridPartFamily;

      struct Traits
      {
        typedef CoordFunction CoordFunctionType;
        typedef typename CoordFunctionType::GridPartType HostGridPartType;

	using LocalCoordFunction = typename CoordFunctionType :: LocalFunctionType;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension - codim, dimensionworld, const GridPartFamily, GeoGeometry > Geometry;
          typedef typename HostGridPartType::template Codim< codim >::LocalGeometryType LocalGeometry;

          typedef Dune::Entity< codim, dimension, const GridPartFamily, GeoEntity > Entity;
          typedef typename HostGridPartType::GridType::template Codim< codim >::EntitySeed EntitySeed;
          typedef Dune::EntityPointer< const GridPartFamily, DefaultEntityPointer< Entity > > EntityPointer;
        };

        typedef DeadIntersection< const GridPartFamily > IntersectionImplType;
        typedef DeadIntersectionIterator< const GridPartFamily > IntersectionIteratorImplType;

        typedef Dune::Intersection< const GridPartFamily, IntersectionImplType > LeafIntersection;
        typedef Dune::Intersection< const GridPartFamily, IntersectionImplType > LevelIntersection;

        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > LeafIntersectionIterator;
        typedef Dune::IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const GridPartFamily, DeadIterator< typename Codim< 0 >::Entity > > HierarchicIterator;
      };

      template< int codim >
      struct Codim
      : public Traits::template Codim< codim >
      {};

      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      typedef typename Traits::HierarchicIterator HierarchicIterator;
    };



    template< class CoordFunction >
    struct GeoGridPartTraits
    {
      typedef GeoGridPart< CoordFunction > GridPartType;
      typedef GeoGridPartFamily< CoordFunction > GridPartFamily;
      typedef GeoGridPartFamily< CoordFunction > GridFamily;

      typedef typename GridPartFamily::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::GridType GridType;

      //! type of twist utility
      typedef MetaTwistUtility< typename HostGridPartType :: TwistUtilityType >  TwistUtilityType;

      typedef IdIndexSet< const GridPartFamily > IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;

      typedef GeoIntersection< const GridPartFamily > IntersectionImplType;
      typedef GeoIntersectionIterator< const GridPartFamily > IntersectionIteratorImplType;

      typedef IntersectionIterator< const GridPartFamily, IntersectionIteratorImplType, IntersectionImplType > IntersectionIteratorType;

      template< int codim >
      struct Codim
      {
        typedef typename GridPartFamily::Traits::template Codim< codim >::Geometry GeometryType;
        typedef typename GridPartFamily::Traits::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridPartFamily::Traits::template Codim< codim >::EntityPointer EntityPointerType;
        typedef typename GridPartFamily::Traits::template Codim< codim >::Entity EntityType;

        typedef typename GridPartFamily::Traits::template Codim< codim >::EntitySeed EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef EntityIterator< codim, const GridPartFamily, GeoIterator< codim, pitype, const GridPartFamily > > IteratorType;
        };
      };

      typedef typename HostGridPartType::CollectiveCommunicationType CollectiveCommunicationType;

      static const bool conforming = HostGridPartType::Traits::conforming;
    };



    // GeoGridPart
    // -----------

    template< class CoordFunction >
    class GeoGridPart
    : public GridPartInterface< GeoGridPartTraits< CoordFunction > >
    {
      typedef GeoGridPart< CoordFunction > ThisType;
      typedef GridPartInterface< GeoGridPartTraits< CoordFunction > > BaseType;

      typedef typename GeoGridPartTraits< CoordFunction >::GridPartFamily GridPartFamily;

      typedef typename GridPartFamily::Traits::HostGridPartType HostGridPartType;

    public:
      typedef CoordFunction CoordFunctionType;

      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::CollectiveCommunicationType CollectiveCommunicationType;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit GeoGridPart ( const CoordFunctionType &coordFunction )
      : coordFunction_( coordFunction ),
        indexSet_( hostGridPart().indexSet() )
      {}

      const GridType &grid () const
      {
        return hostGridPart().grid();
      }

      GridType &grid ()
      {
        return const_cast< GridType & >( hostGridPart().grid() );
      }

      const IndexSetType &indexSet () const
      {
        return indexSet_;
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      begin () const
      {
        return begin< codim, InteriorBorder_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      begin () const
      {
        return GeoIterator< codim, pitype, const GridPartFamily >( coordFunction_, hostGridPart().template begin< codim, pitype >() );
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        return GeoIterator< codim, pitype, const GridPartFamily >( coordFunction_, hostGridPart().template end< codim, pitype >() );
      }

      int level () const
      {
        return hostGridPart().level();
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeoIntersectionIterator< const GridPartFamily >( entity, hostGridPart().ibegin( entity.impl().hostEntity() ) );
      }

      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        return GeoIntersectionIterator< const GridPartFamily >( entity, hostGridPart().iend( entity.impl().hostEntity() ) );
      }

      int boundaryId ( const IntersectionType &intersection ) const
      {
        return hostGridPart().boundaryId( intersection.impl().hostIntersection() );
      }

      const CollectiveCommunicationType &comm () const { return hostGridPart().comm(); }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        GeoDataHandle< GridPartFamily, HostHandleType > handleWrapper( coordFunction_, handle );
        hostGridPart().communicate( handleWrapper, iftype, dir );
      }

      template< class LocalFunction >
      typename Codim< 0 >::EntityType
      exchangeGeometry ( const typename Codim< 0 >::EntityType &entity,
                         const LocalFunction &localCoordFunction ) const
      {
        return typename Codim< 0 >::EntityType::Implementation( entity.impl(), localCoordFunction );
      }

      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return typename Codim< EntitySeed::codimension >::EntityType
                 ::Implementation( coordFunction_, hostGridPart().entity( seed ) );
      }

      // convert a grid entity to a grid part entity ("Gurke!")
      template< class Entity >
      MakeableInterfaceObject< typename Codim< Entity::codimension >::EntityType >
      convert ( const Entity &entity ) const
      {
        // create a grid part entity from a given grid entity
        typedef typename Codim< Entity::codimension >::EntityType EntityType;
        typedef typename EntityType::Implementation Implementation;
        typedef MakeableInterfaceObject< EntityType > EntityObj;

        // here, grid part information can be passed, if necessary
        return EntityObj( Implementation( coordFunction_, entity ) );
      }

      // return reference to the coordfunction
      const CoordFunctionType &coordFunction () const { return coordFunction_; }

    private:
      const HostGridPartType &hostGridPart () const
      {
        return coordFunction_.gridPart();
      }

      const CoordFunctionType &coordFunction_;
      IndexSetType indexSet_;
    };



    // GridEntityAccess for GeoEntity
    // ------------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::Entity< codim, dim, GridFamily, GeoEntity > >
    {
      typedef Dune::Entity< codim, dim, GridFamily, GeoEntity > EntityType;
      typedef GridEntityAccess< typename EntityType::Implementation::HostEntityType > HostAccessType;
      typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return HostAccessType::gridEntity( entity.impl().hostEntity() );
      }
    };



    // EntitySearch for GeoGridPart
    // ----------------------------

    template< class CoordFunction, int codim, PartitionIteratorType partition >
    class EntitySearch< GeoGridPart< CoordFunction >, codim, partition >
    : public DefaultEntitySearch< GeoGridPart< CoordFunction >, codim, partition >
    {
      typedef EntitySearch< GeoGridPart< CoordFunction >, codim, partition > ThisType;
      typedef DefaultEntitySearch< GeoGridPart< CoordFunction >, codim, partition > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit EntitySearch ( const GridPartType &gridPart )
      : BaseType( gridPart )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_HH

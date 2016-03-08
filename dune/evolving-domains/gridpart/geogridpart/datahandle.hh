#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_DATAHANDLE_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_DATAHANDLE_HH

#include <type_traits>

#include <dune/grid/common/datahandleif.hh>

#include <dune/evolving-domains/gridpart/geogridpart/entity.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoDataHandle
    // -------------

    template< class GridFamily, class WrappedHandle >
    class GeoDataHandle
    : public CommDataHandleIF< GeoDataHandle< GridFamily, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      template< class HostEntity >
      class EntityProxy;

    public:
      typedef typename Traits::CoordFunctionType CoordFunctionType;

      GeoDataHandle ( const CoordFunctionType &coordFunction, WrappedHandle &handle )
      : coordFunction_( coordFunction ),
        wrappedHandle_( handle )
      {}

      bool contains ( int dim, int codim ) const
      {
        return wrappedHandle_.contains( dim, codim );
      }

      bool fixedsize ( int dim, int codim ) const
      {
        return wrappedHandle_.fixedsize( dim, codim );
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( coordFunction_, hostEntity );
        return wrappedHandle_.size( *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( coordFunction_, hostEntity );
        wrappedHandle_.gather( buffer, *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        EntityProxy< HostEntity > proxy( coordFunction_, hostEntity );
        wrappedHandle_.scatter( buffer, *proxy, size );
      }

    private:
      const CoordFunctionType &coordFunction_;
      WrappedHandle &wrappedHandle_;
    };



    // GeoDataHandle::EntityProxy
    // --------------------------

    template< class GridFamily, class WrappedHandle >
    template< class HostEntity >
    struct GeoDataHandle< GridFamily, WrappedHandle >::EntityProxy
    {
      static const int dimension = HostEntity::dimension;
      static const int codimension = HostEntity::codimension;

      typedef Dune::Entity< codimension, dimension, const GridFamily, GeoEntity > Entity;

    private:
      typedef GeoEntity< codimension, dimension, const GridFamily > EntityImpl;

    public:
      EntityProxy ( const CoordFunctionType &coordFunction, const HostEntity &hostEntity )
      : entity_( EntityImpl( coordFunction, hostEntity ) )
      {}

      const Entity &operator* () const
      {
        return entity_;
      }

    private:
      Entity entity_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_DATAHANDLE_HH

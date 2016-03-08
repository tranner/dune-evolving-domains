#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH

#include <type_traits>

#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/geometrygrid/geometry.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/evolving-domains/gridpart/geogridpart/cornerstorage.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoGeometryTraits
    // -----------------

    template< class GridFamily >
    class GeoGeometryTraits
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ctype > > MatrixHelper;

      static ctype tolerance () { return 16 * std::numeric_limits< ctype >::epsilon(); }

      template< int mydim, int cdim >
      struct CornerStorage
      {
        typedef GeoCornerStorage< mydim, cdim, GridFamily > Type;
      };

      template< int mydim >
      struct hasSingleGeometryType
      : public Dune::GeoGrid::InferHasSingleGeometryType< GridPartCapabilities::hasSingleGeometryType< HostGridPartType >, dimension, mydim >
      {};
    };



    // GeoGeometry
    // -----------

    template< int mydim, int cdim, class GridFamily >
    class GeoGeometry
    {
      typedef GeoGeometry< mydim, cdim, GridFamily > ThisType;
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;
      typedef typename Traits::CoordFunctionType CoordFunctionType;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;
      using LocalFunctionType = typename CoordFunctionType :: LocalFunctionType;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int codimension = dimension - mydimension;

    private:
      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      typedef typename HostEntityType :: Geometry HostGeometry;

    public:
      //! type of local coordinates
      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      //! type of local coordinates
      typedef FieldVector< ctype, dimension > DimensionCoordinate;
      //! type of global coordinates
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

      //! type of jacobian transposed
      typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;
      //! type of jacobian transposed
      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      //! type of jacobian inverse transposed
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

      using ReferenceElementType = Dune::ReferenceElement< ctype, mydimension >;

      GeoGeometry () = delete;

      GeoGeometry ( const HostGeometry &hostGeo, const LocalFunctionType &mapping,
		    const unsigned int index  = 0 )
	: hostGeo_( hostGeo ), mapping_( mapping ), index_( index ), set_( true )
      {
        assert( int( hostGeo.type().dim() ) == mydimension );
      }

      GeoGeometry ( const ThisType &other )
	: hostGeo_( other.hostGeo_ ), mapping_( other.mapping_ ), index_( other.index_ ),
	  set_( other.set_ )
      {
      }

      const ThisType &operator= ( const ThisType &other )
      {
	std::cerr << "not implemented" << std::endl;
	assert(0);
	return *this;
      }

      void evaluate( const LocalCoordinate& x, GlobalCoordinate& y ) const
      {
	const auto& refEl = Dune::ReferenceElements< ctype, dimension >::general( type() );
	const auto& embedMap = refEl.template geometry< codimension >( index_ );
	DimensionCoordinate xElement = embedMap.global( x );

	mapping_.evaluate( xElement, y );
      }
      void jacobian( const LocalCoordinate& x, Jacobian& mat ) const
      {
	const auto dG = hostGeo_.jacobianTransposed( x );

	const auto& refEl = Dune::ReferenceElements< ctype, dimension >::general( type() );
	const auto& embedMap = refEl.template geometry< codimension >( index_ );
	DimensionCoordinate xElement = embedMap.global( x );

	typename LocalFunctionType :: JacobianRangeType dF;
	mapping_.jacobian( xElement, dF );

	for( unsigned int i = 0; i < coorddimension; ++i )
	  for( unsigned int j = 0; j < mydimension; ++j )
	    {
	      mat[ i ][ j ] = 0;
	      for( unsigned int k = 0; k < coorddimension; ++k )
		mat[ i ][ j ] += dF[ i ][ k ] * dG[ j ][ k ];
	    }
      }

      /*operator bool () const { return set_; }*/

      bool set() const { return set_; }

      bool affine () const
      {
	return mapping_.order();
      }
      GeometryType type () const
      {
	return hostGeo_.type();
      }

      int corners () const
      {
	return referenceElement().size( mydimension );
      }
      GlobalCoordinate corner ( const int i ) const
      {
	GlobalCoordinate y;
	evaluate( referenceElement().position( i, mydimension ), y );
	return y;
      }
      GlobalCoordinate center () const
      {
	std::cerr << "warning! not implemented!" << std::endl;
	GlobalCoordinate y;
	// mapping_.evaluate( referenceElement().position( 0, 0 ), y );
	return y;
      }

      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
	GlobalCoordinate y;
	mapping_.evaluate( local, y );
	return y;
      }
      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
	std::cerr << "not implemented" << std::endl;
	assert(0);
	return LocalCoordinate();
      }

      ctype integrationElement ( const LocalCoordinate &local ) const
      {
	Jacobian jac;
	jacobian( local, jac );

	Dune::FieldMatrix< ctype, mydimension, mydimension > jjt;
	for( unsigned int i = 0; i < mydimension; ++i )
	  for( unsigned int j = 0; j < mydimension; ++j )
	    {
	      jjt[ i ][ j ] = 0;
	      for( unsigned int k = 0; k < coorddimension; ++k )
		jjt[ i ][ j ] += jac[ k ][ i ] * jac[ k ][ j ];
	    }

	const double det = jjt.determinant();
	assert( det > 0 );

	return std::sqrt( det );
      }
      ctype volume () const
      {
	std::cerr << "not implemented" << std::endl;
	assert(0);
	return -1;
      }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
      {
	Jacobian jac;
	jacobian( local, jac );

	JacobianTransposed jact;
	for( unsigned int i = 0; i < mydimension; ++i )
	  for( unsigned int j = 0; j < coorddimension; ++j )
	    jact[ i ][ j ] = jac[ j ][ i ];

	return jact;
      }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
	Jacobian jac;
	jacobian( local, jac );

	Dune::FieldMatrix< ctype, mydimension, mydimension > jjt;
	for( unsigned int i = 0; i < mydimension; ++i )
	  for( unsigned int j = 0; j < mydimension; ++j )
	    {
	      jjt[ i ][ j ] = 0;
	      for( unsigned int k = 0; k < coorddimension; ++k )
		jjt[ i ][ j ] += jac[ k ][ i ] * jac[ k ][ j ];
	    }

	jjt.invert();

	JacobianInverseTransposed ret;
	for( unsigned int i = 0; i < mydimension; ++i )
	  for( unsigned int j = 0; j < coorddimension; ++j )
	    {
	      ret[ j ][ i ] = 0.0;
	      for( unsigned int k = 0; k < mydimension; ++k )
		ret[ j ][ i ] += jjt[ i ][ k ] * jac[ j ][ k ];
	    }

	return ret;
      }

    protected:
      const ReferenceElementType& referenceElement() const
      {
	return Dune::ReferenceElements< ctype, mydimension >::general( type() );
      }

    private:
      const HostGeometry hostGeo_;
      const LocalFunctionType mapping_;
      const unsigned int index_;
      bool set_;
    };
  } // namespace Fem


} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH

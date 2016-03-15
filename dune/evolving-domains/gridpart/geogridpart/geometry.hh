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
      using LocalFunctionType = typename CoordFunctionType :: LocalFunctionType;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int codimension = dimension - mydimension;

    private:
      typedef typename Traits::HostGridPartType HostGridPartType;
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      typedef typename HostEntityType :: Geometry HostGeometry;

      typedef typename HostGridPartType::template Codim< 0 >::EntityType HostEntity0Type;
      typedef typename HostEntity0Type :: EntitySeed HostEntity0SeedType;


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

      GeoGeometry( const CoordFunctionType& coordFunction, const HostEntity0SeedType& seed,
		   const HostGeometry& hostGeometry, const unsigned int index = 0 )
	: coordFunction_( &coordFunction ), seed_( seed ),
	  hostGeometry_( hostGeometry ), index_( index )
      {}

      GeoGeometry( const ThisType& other )
	: coordFunction_( other.coordFunction_ ), seed_( other.seed_ ),
	  hostGeometry_( other.hostGeometry_ ), index_( other.index_ )
      {}

      const ThisType &operator= ( const ThisType &other )
      {
	std::cerr << "not implemented" << std::endl;
	assert(0);
	return *this;
      }

      void evaluate( const LocalCoordinate& x, GlobalCoordinate& y ) const
      {
	Dune::GeometryType geoType;
	geoType.makeSimplex( dimension );

	const auto& refEl = Dune::ReferenceElements< ctype, dimension >::general( geoType );
	const auto& embedMap = refEl.template geometry< codimension >( index_ );
	DimensionCoordinate xElement = embedMap.global( x );

	HostEntity0Type e = coordFunction().gridPart().grid().entity( seed_ );
	const LocalFunctionType lf = coordFunction().localFunction( e );
	lf.evaluate( xElement, y );
      }
      void jacobian( const LocalCoordinate& x, Jacobian& mat ) const
      {
	const auto dG = hostGeometry_.jacobianTransposed( x );

	Dune::GeometryType geoType;
	geoType.makeSimplex( dimension );

	const auto& refEl = Dune::ReferenceElements< ctype, dimension >::general( geoType );
	const auto& embedMap = refEl.template geometry< codimension >( index_ );
	DimensionCoordinate xElement = embedMap.global( x );

	typename LocalFunctionType :: JacobianRangeType dF;
	HostEntity0Type e = coordFunction().gridPart().grid().entity( seed_ );
	const LocalFunctionType lf = coordFunction().localFunction( e );

	lf.jacobian( xElement, dF );

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
	const HostEntity0Type e = coordFunction().gridPart().grid().entity( seed_ );
	const LocalFunctionType lf = coordFunction().localFunction( e );
	return lf.order();
      }
      GeometryType type () const
      {
	return hostGeometry_.type();
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
	const HostEntity0Type e = coordFunction().gridPart().grid().entity( seed_ );
	const LocalFunctionType lf = coordFunction().localFunction( e );
	lf.evaluate( local, y );
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
      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

      const ReferenceElementType& referenceElement() const
      {
	return Dune::ReferenceElements< ctype, mydimension >::general( type() );
      }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
      const HostEntity0SeedType seed_;
      const HostGeometry hostGeometry_;
      const unsigned int index_;
      bool set_;
    };
  } // namespace Fem


} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_GEOMETRY_HH

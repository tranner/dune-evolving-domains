#ifndef ELLIPTIC_COUPLED_HH
#define ELLIPTIC_COUPLED_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/grid/common/entitypointer.hh>

// EllipticOperator
// ----------------

//! [Class for elliptic operator]
template< class DiscreteFunction, class Model, unsigned int cd >
struct EllipticOperator
  : public virtual Dune::Fem::Operator< DiscreteFunction >
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model ModelType;

protected:

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionType::DomainType GlobalCoordType;

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType :: Intersection IntersectionType;
  typedef typename IntersectionType :: Geometry IntersectionGeometryType;

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef typename QuadratureType::CoordinateType LocalCoordType;

public:
  EllipticOperator ( const ModelType &model )
    : model_( model ), codim( cd )
  {}

  // prepare the solution vector
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u )
  {
    // nothing to do -- closed surfaces only
  }


  //! application operator
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

protected:
  const ModelType &model() const { return model_; }

private:
  ModelType model_;

protected:
  const unsigned int codim;
};

// DifferentiableEllipticOperator
// ------------------------------

//! [Class for linearizable elliptic operator]
template< class JacobianOperator, class Model, unsigned int cd >
struct DifferentiableEllipticOperator
  : public EllipticOperator< typename JacobianOperator :: DomainFunctionType,
			     Model, cd >,
    public Dune :: Fem :: DifferentiableOperator< JacobianOperator >
{
  typedef EllipticOperator< typename JacobianOperator::DomainFunctionType, Model, cd > BaseType;

  typedef JacobianOperator JacobianOperatorType;

  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;

protected:
  typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename BaseType::LocalFunctionType LocalFunctionType;

  typedef typename BaseType::IteratorType IteratorType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::GeometryType GeometryType;

  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename BaseType::IntersectionType IntersectionType;
  typedef typename BaseType::IntersectionGeometryType IntersectionGeometryType;

  typedef typename BaseType::GlobalCoordType GlobalCoordType;
  typedef typename BaseType::LocalCoordType LocalCoordType;

  typedef typename BaseType::QuadratureType QuadratureType;
  typedef typename BaseType::FaceQuadratureType FaceQuadratureType;

public:
  //! constructor
  DifferentiableEllipticOperator( const ModelType &model )
    : BaseType( model )
  {}

  //! method to setup the jacobian of the operator for storage in a matrix
  void jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const;

protected:
  using BaseType::model;
  using BaseType::codim;
};

// Implementation of EllipticOperator
// ----------------------------------

template< class DiscreteFunction, class Model, unsigned int cd >
void EllipticOperator< DiscreteFunction, Model, cd >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
{
  // clear destination
  w.clear();

  // get discrete function space
  const DiscreteFunctionSpaceType &dfSpace = w.space();

  // iterate over grid
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      // get entity
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      // get local representation of the discrete functions
      const LocalFunctionType uLocal = u.localFunction( entity );
      LocalFunctionType wLocal = w.localFunction( entity );

      // obtain quadrature order
      const int quadOrder = uLocal.order() + wLocal.order();

      QuadratureType quadrature( entity, quadOrder );
      size_t numQuadraturePoints = quadrature.nop();

      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  // obatain quadrature points in quadrature coords
	  const LocalCoordType &xLocal = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) *
	    geometry.integrationElement( xLocal );

	  // evaluate discrete solutions and mass flux
	  typename LocalFunctionType::RangeType uhx, cuhx;
	  uLocal.evaluate( xLocal, uhx );
	  model().source( entity, xLocal, uhx, cuhx );

	  // evaluate discrete gradient and diffusive flux
	  typename LocalFunctionType::JacobianRangeType duhx, aduhx;
	  uLocal.jacobian( xLocal, duhx );
	  model().diffusiveFlux( entity, xLocal, cuhx, duhx, aduhx );

	  if( cd == 1 )
	    {
	      cuhx += 1.0 * uhx;
#warning using beta = 1 here
	    }

	  // multiply by quadrature weight
	  cuhx *= weight;
	  aduhx *= weight;

	  // add cuhx * phi and aduhx \cdot grad phi_i to wLocal[ i ]
	  wLocal.axpy( quadrature[ pt ], cuhx, aduhx );
	}

      if( cd == 0 )
	{
	  const IntersectionIteratorType iEnd = dfSpace.gridPart().iend( entity );
	  for( IntersectionIteratorType iIt = dfSpace.gridPart().ibegin( entity );
	       iIt != iEnd; ++ iIt )
	    {
	      if( iIt->boundary() )
		{
		  // get intersection
		  const IntersectionType &intersection = *iIt;
		  const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

		  FaceQuadratureType faceQuadrature( dfSpace.gridPart(), intersection,
						     quadOrder, FaceQuadratureType::INSIDE );
		  size_t numFaceQuadraturePoints = faceQuadrature.nop();

		  for( size_t pt = 0; pt < numFaceQuadraturePoints; ++pt )
		    {
		      // obtain quadrature point
		      const typename FaceQuadratureType :: LocalCoordinateType &xQuad = faceQuadrature.localPoint( pt );
		      const double weight = faceQuadrature.weight( pt ) * intersectionGeometry.integrationElement( xQuad );

		      // evaluate discrete function and flux
		      typename LocalFunctionType :: RangeType uhx, cuhx;
		      uLocal.evaluate( faceQuadrature.point( pt ), uhx );
		      cuhx = 1.0 * uhx;
#warning using alpha = 1 here

		      // multiply by quadrature weight
		      cuhx *= weight;

		      // add cuhx * phi to wLocal[ i ]
		      wLocal.axpy( faceQuadrature[ pt ], cuhx );
		    }
		}
	    }
	}
    }

  // communicate in parallel runs
  w.communicate();
}

// Implementation of DifferentiableEllipticOperator
// ------------------------------------------------

template < class JacobianOperator, class Model, unsigned int cd >
void DifferentiableEllipticOperator< JacobianOperator, Model, cd >
::jacobian ( const DiscreteFunctionType &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  // get discrete function space
  const DiscreteFunctionSpaceType &dfSpace = u.space();

  // set up matrix stencil
  Dune :: Fem :: DiagonalStencil< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > stencil( dfSpace, dfSpace );
  jOp.reserve( stencil );
  jOp.clear();

  // set up basis function storage
  const int blockSize = dfSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename LocalFunctionType::RangeType > phi( dfSpace.blockMapper().maxNumDofs()*blockSize );
  std::vector< typename LocalFunctionType::JacobianRangeType > dphi( dfSpace.blockMapper().maxNumDofs()*blockSize );

  // loop over grid
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      // find entity
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      // construct local representations of data
      const LocalFunctionType uLocal = u.localFunction( entity );
      LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

      // find local basis functions
      const BasisFunctionSetType &basisSet = jLocal.domainBasisFunctionSet();
      const unsigned int numBasisFunctions = basisSet.size();

      // perform quadrature loop
      QuadratureType quadrature( entity, 2*dfSpace.order() );
      size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  // obatain quadrature points in quadrature coords
	  const LocalCoordType &xLocal = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) *
	    geometry.integrationElement( xLocal );

	  // evaluate basis functions at quadrature points
	  basisSet.evaluateAll( xLocal, phi );
	  basisSet.jacobianAll( xLocal, dphi );

	  // get value for linearization
	  typename LocalFunctionType :: RangeType ux;
	  typename LocalFunctionType :: JacobianRangeType dux;
	  uLocal.evaluate( quadrature[ pt ], ux );
	  uLocal.jacobian( quadrature[ pt ], dux );

	  for( unsigned int i = 0; i < numBasisFunctions; ++i )
	    {
	      // evaluate fluxes
	      typename LocalFunctionType::RangeType cphi;
	      model().linSource( ux, entity, xLocal, phi[ i ], cphi );

	      if( cd == 1 )
		{
		  cphi += 1.0 * phi[ i ];
#warning using beta = 1 here
		}

	      typename LocalFunctionType::JacobianRangeType adphi;
	      model().linDiffusiveFlux( ux, dux, entity, xLocal, phi[ i ], dphi[ i ], adphi );

	      // add contribution to local matrix
	      jLocal.column( i ).axpy( phi, dphi, cphi, adphi, weight );
	    }
	}

      if( cd == 0 )
	{ // then coupling terms if bulk surface terms
	  const IntersectionIteratorType iEnd = dfSpace.gridPart().iend( entity );
	  for( IntersectionIteratorType iIt = dfSpace.gridPart().ibegin( entity );
	       iIt != iEnd; ++iIt )
	    {
	      if( iIt->boundary() )
		{
		  // get intersection
		  const IntersectionType &intersection = *iIt;
		  const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

		  FaceQuadratureType faceQuadrature( dfSpace.gridPart(), intersection,
						     2.0*dfSpace.order(),
						     FaceQuadratureType::INSIDE );
		  size_t numFaceQuadraturePoints = faceQuadrature.nop();

		  for( size_t pt = 0; pt < numFaceQuadraturePoints; ++pt )
		    {
		      // obtain quadrature point
		      const typename FaceQuadratureType :: LocalCoordinateType &xQuad = faceQuadrature.localPoint( pt );
		      const double weight = faceQuadrature.weight( pt ) * intersectionGeometry.integrationElement( xQuad );

		      // evaluate basis functions at quadrature  points
		      basisSet.evaluateAll( faceQuadrature.point( pt ), phi );

		      for( unsigned int i = 0; i < numBasisFunctions; ++i )
			{
			  // evaluate coupling flux
			  typename LocalFunctionType :: RangeType cphi;
			  cphi = 1.0 * phi[ i ];
#warning using alpha = 1 here

			  // add contribution to local matrix
			  jLocal.column( i ).axpy( phi, cphi, weight );
			}
		    }
		}
	    }
	}
    }
  jOp.communicate();
}

#endif // #ifndef ELLIPTIC_COUPLED_HH

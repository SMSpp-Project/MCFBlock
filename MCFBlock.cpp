/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MCFBlock class.
 *
 * \version 0.30
 *
 * \date 21 - 03 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// in the DIMACS format, comment lines start with 'c'
static inline std::istream & eatDMXcomments( std::istream& is )
{
 for(;;) {
  is >> std::ws;  // skip whitespaces
  if( is.peek() == is.widen( 'c' ) )
   // a comment: skip the rest of line and move to next
   is.ignore( std::numeric_limits<std::streamsize>::max() , is.widen( '\n' )
	      );
  else
   break;
  }

 return( is );
 }

/*--------------------------------------------------------------------------*/

static inline void print_UB( std::ostream& os , MCFBlock::FNumber ub )
{
 if( ub == Inf<MCFBlock::FNumber>() )
  os << "+Inf";
 else
  os << ub;
 }

/*--------------------------------------------------------------------------*/

static MCFBlock::FNumber read_UB( std::istream & iStrm )
{
 iStrm >> eatcomments;
 int c = iStrm.peek();
 if( ! iStrm )
  throw( std::invalid_argument( "error reading the input stream" ) );
  
 if( ( c != 'I' ) && ( c != 'i' ) ) {
  MCFBlock::FNumber res;
  iStrm >> res;
  if( ! iStrm )
   throw( std::invalid_argument( "error reading the input stream" ) );
  return( res );
  }

 do { c = iStrm.get(); c = iStrm.peek();
      if( ! iStrm )
       throw( std::invalid_argument( "error reading the input stream" ) );

  } while( ( c != iStrm.widen( ' ' ) ) &&
	   ( c != iStrm.widen( '\n' ) ) &&
	   ( c != iStrm.widen( '\t' ) ) );

 return( Inf<MCFBlock::FNumber>() );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ
template< typename T >
static inline MCFBlock::Index countdiff( std::vector<T>:const_iterator beg ,
					 std::vector<T>:const_iterator end ,
					 std::vector<T>:const_iterator cmp )
{
 MCFBlock::Index ndiff = 0;
 for( ; beg != end ; )
  if( *(beg++) != *(cmp++) )
   ndiff++;

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices
template< typename T >
static inline MCFBlock::Index countdiff( std::vector<T> & vec ,
					 MCFBlock::c_Vec_Index & nms ,
					 std::vector<T>:const_iterator cmp ,
					 MCFBlock::c_Index n_max )
{
 MCFBlock::Index ndiff = 0;
 for( auto beg = nms.begin() ; beg != nms.end() ; ) {
  if( nms[ *beg ] >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ *(beg++) ] != *(cmp++) )
   ndiff++;
  }

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another
template< typename T >
static inline MCFBlock::Index copyidx( std::vector<T> & vec ,
				       MCFBlock::c_Vec_Index & nms ,
				       std::vector<T>:const_iterator cpy )
{
 for( auto nit = nms.begin() ; nit < nms.end() ; )
  vec[ *(nit++) ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register MCFBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( MCFBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF MCFBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::load( c_Index n , c_Vec_Index & pEn , c_Vec_Index & pSn ,
		     c_Vec_FNumber & pU , c_Vec_CNumber & pC ,
		     c_Vec_FNumber & pB , c_Index dn , c_Index dm ,
		     c_Index mdn , c_Index mdm )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( pEn.size() != pSn.size() )
  throw( std::invalid_argument( "pEn size != pSn size" ) );

 if( ( pC.size() > 0 ) && ( pC.size() != pSn.size() ) )
  throw( std::invalid_argument( "pCn size != pSn size" ) );

 if( ( pU.size() > 0 ) && ( pU.size() != pSn.size() ) )
  throw( std::invalid_argument( "pUn size != pSn size" ) );

 if( ( pB.size() > 0 ) && ( pB.size() != n ) )
  throw( std::invalid_argument( "pB size != n" ) );

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( NNodes )
  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNodes = n;
 NArcs = pSn.size();
 MaxNNodes = NNodes + ( mdn > dn ? mdn - dn : 0 );
 c_Index MaxNArcs = NArcs + ( mdm > dm ? mdm - dm : 0 );
 NStaticNodes = dn > n : 0 ? n - dn;
 NStaticArcs = dm > m : 0 ? m - dm;

 SN.resize( MaxNArcs , 0 );
 std::copy( pSn.begin() , pSn.end() , SN.begin() );

 EN.resize( get_MaxNArcs() , 0 );
 std::copy( pEn.begin() , pEn.end() , EN.begin() );

 if( std::any_of( pC.begin() , pC.end() ,
		  []( c_CNumber ci ) { return( ci != 0 ); } ) ) {
  C.resize( get_MaxNArcs() , 0 );
  std::copy( pC.begin() , pC.end() , C.begin() );
  }

 if( std::any_of( pU.begin() , pU.end() ,
		  []( c_FNumber ui ) { return( ui < Inf<FNumber>() ); } ) ) {
  U.resize( get_MaxNArcs() , Inf<FNumber>() );
  std::copy( pU.begin() , pU.end() , U.begin() );
  }

 if( std::any_of( pB.begin() , pB.end() ,
		  []( c_FNumber bi ) { return( bi != 0 ); } ) ) {
  B.resize( get_MaxNNodes() , 0 );
  std::copy( pB.begin() , pB.end() , B.begin() );
  }

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 }  // end( MCFBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::load( std::istream &input )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( NNodes )
  guts_of_destructor();

 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -

 char c;
 if( ! ( input >> eatDMXcomments >> c ) )
  throw( std::invalid_argument( "error reading the input stream" ) );

 if( c != 'p' )
  throw( std::invalid_argument( "format error in the input stream" ) );

 input >> eatDMXcomments;
 input.ignore( 3 , ' ' );  // skip "min"

 if( ! ( input >> eatDMXcomments >> NNodes ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of nodes" ) );

 Index tm;
 if( ! ( input >> eatDMXcomments >> NArcs ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of arcs" ) );

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SN.resize( NArcs );
 EN.resize( NArcs );
 C.resize( NArcs , 0 );
 U.resize( NArcs , Inf<FNumber>() );
 B.resize( NNodes , 0 );

 NStaticNodes = MaxNNodes = NNodes;
 NStaticArcs = NArcs;

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i = 0;  // arc counter
 for(;;) {
  if( ! ( input >> eatDMXcomments >> c ) )  // read next descriptor
   break;                                   // if none, end

  switch( c ) {
   case( 'n' ):  // description of a node
    Index j;
    if( ! ( input >> j ) )
     throw( std::invalid_argument( "error reading node name" ) );

    if( ( j < 1 ) || ( j > NNodes ) )
     throw( std::invalid_argument( "invalid node name" ) );

    FNumber Dfctj;
    if( ! ( input >> Dfctj ) )
     throw( std::invalid_argument( "error reading deficit" ) );

    B[ j - 1 ] -= Dfctj;
    break;

   case( 'a' ):  // description of an arc
    if( i == NArcs )
     throw( std::invalid_argument( "too many arc descriptors" ) );

    if( ! ( input >> SN[ i ] ) )
     throw( std::invalid_argument( "error reading start node" ) );

    if( ( SN[ i ] < 1 ) || ( SN[ i ] > NNodes ) )
     throw( std::invalid_argument( "invalid start node" ) );

    if( ! ( input >> EN[ i ] ) )
     throw( std::invalid_argument( "error reading end node" ) );

    if( ( EN[ i ] < 1 ) || ( EN[ i ] > NNodes ) )
     throw( std::invalid_argument( "LoadDMX: invalid end node" ) );

    if( SN[ i ] == EN[ i ] )
     throw( std::invalid_argument( "self-loops not permitted" ) );

    FNumber LB;
    if( ! ( input >> LB ) )
     throw( std::invalid_argument( "error reading lower bound" ) );

    U[ i ] = read_UB( input );

    if( ! ( input >> C[ i ] ) )
     throw( std::invalid_argument( "error reading arc cost" ) );

    if( U[ i ] < LB )
     throw( std::invalid_argument( "lower bound > upper bound" ) );

    if( LB > 0 ) {
     if( U[ i ] < Inf<MCFBlock::FNumber>() )
      U[ i ] -= LB;
     B[ SN[ i ] - 1 ] += LB;
     B[ EN[ i ] - 1 ] -= LB;
     }
    i++;
    break; 

   default:  // invalid code- - - - - - - - - - - - - - - - - - - - - - - - -
    throw( std::invalid_argument( "invalid DMX code" ) );

   }  // end( switch( c ) )
  }  // end( for( ever ) )

 if( i < NArcs )
  throw( std::invalid_argument( "too few arc descriptors" ) );

 // simplify out the deta structures- - - - - - - - - - - - - - - - - - - - -

 if( std::all_of( C.begin() , C.end() ,
		  []( c_CNumber ci ) { return( ci == 0 ); } ) )
  C.clear();

 if( std::all_of( B.begin() , B.end() ,
		  []( c_FNumber bi ) { return( bi == 0 ); } ) )
  B.clear();

 if( std::all_of( U.begin() , U.end() ,
		  []( c_FNumber ui ) { return( ui == Inf<FNumber>() ); } ) )
  U.clear();
 
 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 }  // end( MCFBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::deserialize( netCDF::NcGroup && group , Block * father )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( NNodes )
  guts_of_destructor();
		   
 // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.getDim( "NNodes" );
 if( nn.isNull() )
  throw( std::logic_error( "NNodes dimension is required" ) );
 NNodes = nn.getSize();

 netCDF::NcDim na = group.getDim( "NArcs" );
 if( na.isNull() )
  throw( std::logic_error( "NArcs dimension is required" ) );
 NArcs = na.getSize();

 Index DynNNodes = 0;
 netCDF::NcDim dn = group.getDim( "DynNNodes" );
 if( dn.isNull() )
  NStaticNodes = NNodes;
 else {
  DynNNodes = dn.getSize();
  NStaticNodes = DynNNodes > NNodes ? 0 : NNodes - DynNNodes;
  }

 Index DynNArcs = 0;
 netCDF::NcDim dm = group.getDim( "DynNArcs" );
 if( dm.isNull() )
  NStaticArcs = NArcs;
 else {
  DynNArcs = dm.getSize();
  NStaticArcs = DynNArcs > NArcs ? 0 : NArcs - DynNArcs;
  }

 MaxNNodes = NNodes;
 netCDF::NcDim mdn = group.getDim( "MaxDynNNodes" );
 if( ( ! mdn.isNull() ) && ( mdn.getSize() > DynNNodes ) )
  MaxNNodes += mdn.getSize() - DynNNodes;

 Index MaxNArcs = NArcs;
 netCDF::NcDim mdm = group.getDim( "MaxDynNArcs" );
 if( ( ! mdm.isNull() ) && ( mdm.getSize() > DynNArcs ) )
  MaxNNodes += mdm.getSize() - DynNArcs;
 
 netCDF::NcVar sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::logic_error( "Starting Nodes not found" ) );

 SN.resize( MaxNArcs );

 std::vector<size_t> start = { 0 };
 std::vector<size_t> counta = { NArcs };
 sn.getVar( start , counta , SN.data() );

 netCDF::NcVar en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::logic_error( "Ending Nodes not found" ) );

 EN.resize( MaxNArcs );
 en.getVar( start , counta , EN.data() );

 netCDF::NcVar cst = group.getVar( "C" );
 if( ! cst.isNull() ) {
  C.resize( MaxNArcs , 0 );
  cst.getVar( start , counta , C.data() );
  if( std::all_of( C.begin() , C.begin() + NArcs ,
		   []( c_CNumber ci ) { return( ci == 0 ); } ) )
   C.clear();
  }

 netCDF::NcVar cap = group.getVar( "U" );
 if( ! cap.isNull() ) {
  U.resize( MaxNArcs , Inf<FNumber>() );
  cap.getVar( start , counta , U.data() );
  if( std::all_of( U.begin() , U.begin() + NArcs ,
		   []( c_FNumber ui ) { return( ui == Inf<FNumber>() ); } ) )
   U.clear();
  }

 netCDF::NcVar dfc = group.getVar( "B" );
 if( ! dfc.isNull() ) {
  B.resize( MaxNNodes , 0 );
  std::vector<size_t> countn = { NNodes };
  dfc.getVar( start , countn , B.data() );
  if( std::all_of( B.begin() , B.begin() + NNodes ,
		   []( c_FNumber bi ) { return( bi == 0 ); } ) )
   B.clear();
  }

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 }  // end( MCFBlock::deserialize )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_abstract_variables( Configuration *stvv )
{
 if( x.size() || dx.size() )  // the variables are there already
  return;                     // nothing to do

 if( HasStaticX() ) {
  x.resize( get_NStaticArcs() );
  for( auto & var : x ) {
   var.is_positive( true , eNoBlck );
   var.set_Block( this );
   }

  add_static_variable( x );
  }

 if( MayHaveDynX() ) {
  dx.resize( get_NArcs() - get_NStaticArcs() );
  for( auto & var : dx ) {
   var.is_positive( true , eNoBlck );
   var.set_Block( this );
   }

  add_dynamic_variable( dx );
  } 
 }  // end( MCFBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_abstract_constraints( Configuration *stcc )
{
 if( E.size() || dE.size() )  // the constraints are there already
  return;                     // nothing to do

 // count number of nonzeroes in each constraint, i.e., #FS( i ) + #BS( i )
 std::vector<Index> count( get_NNodes() );
 
 for( Index i = 0 ; i < get_NArcs() ; ++i ) {
  count[ SN[ i ] - 1 ]++;
  count[ EN[ i ] - 1 ]++;
  }

 // initialize the vectors of coefficients, and reset count[]
 std::vector< LinearFunction::v_coeff_pair > coeffs( get_NNodes() );

 for( Index i = 0 ; i < get_NNodes() ; ++i ) {
  coeffs[ i ].resize( count[ i ] );
  count[ i ] = 0;
  }

 // construct the vector of coefficients, static phase
 Index i = 0;
 for( ; i < get_NStaticArcs() ; ++i ) {
  coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( -1 ) );
  coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( 1 ) );
  }

 // construct the vector of coefficients, dynamic phase
 if( MayHaveDynX() )
  for( auto dxi = dx.begin() ; i < get_NArcs() ; ++i , ++dxi ) {
   coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &(*dxi) , double( -1 ) );
   coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &(*dxi) , double( 1 ) );
   }

 // generate the node-arc incidence matrix- - - - - - - - - - - - - - - - - -
 // each constraint is an equality, i.e., LHS = RHS = B[ i ]

 // static part
 if( HasStaticE() ) {
  E.resize( get_NStaticNodes() );

  for( Index i = 0 ; i < get_NStaticNodes() ; ++i ) {
   E[ i ].set_both( B.size() ? B[ i ] : 0 );
   // note that the pairs are not ordered if there is a dynamic part
   E[ i ].set_function( new LinearFunction( std::move( coeffs[ i ] ) , 0 ,
					    ! HasDynamicX() ) );
   E[ i ].set_Block( this );  // this is done last ==> no Modification
   }

  add_static_constraint( E );
  }

 // dynamic part
 if( MayHaveDynE() ) {
  dE.resize( get_NNodes() - get_NStaticNodes() );

  Index i = get_NStaticNodes();
  for( auto & cnst : dE ) {
   cnst.set_both( B.size() ? B[ i ] : 0 );
   // note that the pairs are not ordered since there is a dynamic part
   cnst.set_function( new LinearFunction( std::move( coeffs[ i++ ] ) , 0 ,
					  false ) );
   cnst.set_Block( this );  // this is done last ==> no Modification
   }

  add_dynamic_constraint( dE );
  }

 // generate the bound constraints- - - - - - - - - - - - - - - - - - - - - -
 // if upper bounds are not there, the LB0Constraint can be skipped entirely
 // if the Configuration agrees

 if( U.empty() ) {
  auto tstcc = dynamic_cast<SimpleConfiguration<int> *>( stcc );

  if( ( ! tstcc ) && f_BlockConfig &&
      f_BlockConfig->f_static_constraints_Configuration )
   tstcc = dynamic_cast<SimpleConfiguration<int> *>(
			 f_BlockConfig->f_static_constraints_Configuration );
  if( tstcc && ( tstcc->f_value != 0 ) )
   return;
  }

 // static part
 if( HasStaticX() ) {
  UB.resize( get_NStaticArcs() );
  for( Index i = 0 ; i < get_NStaticArcs() ; ++i ) {
   UB[ i ].set_variable( & x[ i ] , eNoBlck );
   UB[ i ].set_rhs( U[ i ] , eNoBlck );
   UB[ i ].set_Block( this );  // this is done last ==> no Modification
   }

  add_static_constraint( &LB );
  }

 // dynamic part
 if( MayHaveDynX() ) {
  dUB.resize( get_NArcs() - get_NStaticArcs() );

  auto dxi = dx.begin();
  auto ui = U.begin() + get_NStaticArcs();
  for( auto & cnst : dUB ) {
   cnst.set_variable( &(*(dxi++)) , eNoBlck );
   cnst.set_rhs( *(ui++) , eNoBlck );
   cnst.set_Block( this );  // this is done last ==> no Modification
   }

  add_dynamic_constraint( &dUB );
  }
 }  // end( MCFBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_objective( Configuration *objc )
{
 if( ! get_objective().empty() )  // an objective is there already
  return;                         // cowardly (and silently) return

 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 LinearFunction::v_coeff_pair p;
 /*!!
 Index nzc = 0;  // counter of nonzero coefficients
 for( Index i = 0 ; i < C.size() ; i++ )
  if( C[ i ] != 0 )
   nzc++;

 double sprs = 0;    // use NNConstraint rather than BoxConstraint

 auto tobjc = dynamic_cast<SimpleConfiguration<double> *>( objc );

 if( ( ! tobjc ) && f_BlockConfig
                 && f_BlockConfig->f_objective_Configuration )
  tobjc = dynamic_cast<SimpleConfiguration<double> *>(
			          f_BlockConfig->f_objective_Configuration );
 if( tobjc )
  sprs = std::max( Function::FunctionValue( 0 ) ,
		 std::min( tobjc->f_value , Function::FunctionValue( 1 ) ) );

 if( nzc >= std::ceil( sprs * get_NArcs() ) ) !!*/ {
  // construct a "dense" LinearFunction - - - - - - - - - - - - - - - - - - -
  p.resize( get_NArcs() );
  Index i = 0;

  // static part
  if( HasStaticX() )
   if( C.size() )
    for( ; i < get_NStaticArcs() ; ++i ) {
     p[ i ].first = &x[ i ];
     p[ i ].second = C[ i ];
     }
   else
    for( ; i < get_NStaticArcs() ; ++i ) {
     p[ i ].first = &x[ i ];
     p[ i ].second = 0;
     }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   if( C.size() )
    for( ; i < get_NArcs() ; ++i ) {
     p[ i ].first = &(*(dxi++));
     p[ i ].second = C[ i ];
     }
   else
    for( ; i < get_NArcs() ; ++i ) {
     p[ i ].first = &(*(dxi++));
     p[ i ].second = 0;
     }
   }
  }
 /*!!
 else {
  // construct a "sparse" LinearFunction- - - - - - - - - - - - - - - - - - -
  p.resize( nzc );
  for( Index i = 0 , j = 0 ; i < C.size() ; ++i )
   if( C[ i ] ) {
    p[ j ].first = &x[ i ];
    p[ j++ ].second = C[ i ];
    }
  }
  !!*/

 // ensure no Modification is issued: this may happen in case a MCFBlock
 // is re-loaded, so that set_objective( c ) had already been called
 // note that the pairs are not ordered if there is a dynamic part
 c.set_function( new LinearFunction( std::move( p ) , 0 , ! HasDynamicX() ) ,
		 eNoMod );
 c.set_Block( this );

 set_objective( c , eNoMod );

 }  // end( MCFBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool MCFBlock::flow_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  // static part
  for( const auto & cnst : E )
   if( cnst.rel_viol() > feps )
    return( false );

  // dynamic part
  for( const auto & cnst : dE )
   if( cnst.rel_viol() > feps )
    return( false );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Index i = 0;
  Vec_FNumber tB = B;

  // static part
  for(  ; i < get_NStaticArcs() ; ++i ) {
   c_FNumber xi = x[ i ].get_value();
   tB[ SN[ i ] - 1 ] += xi;
   tB[ EN[ i ] - 1 ] -= xi;
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for(  ; i < get_NArcs() ; ++i ) {
    c_FNumber xi = (*(dxi++)).get_value();
    tB[ SN[ i ] - 1 ] += xi;
    tB[ EN[ i ] - 1 ] -= xi;
    }
   }

 for( Index i = 0 ; i < get_NNodes() ; ++i ) {
   c_FNumber slck = B[ i ] == 0 ? std::abs( tB[ i ] ) :
                                  std::abs( tB[ i ] / B[ i ] );
   if( slck > feps )
    return( false );
   }
  }

 return( true );

 }  // end( MCFBlock::flow_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::bound_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  // static part
  if( HasStaticX() )
   if( UB.empty() ) {
    for( const auto & var : x )
     if( ! var.is_feasible() )
      return( false );
    }
   else
    for( const auto & cnst : UB )
     if( cnst.rel_viol() > feps )
      return( false );

  // dynamic part
  if( HasDynamicX() )
   if( dUB.empty() ) {
    for( const auto & var : dx )
     if( ! var.is_feasible() )
      return( false );
    }
   else
    for( const auto & cnst : dUB )
     if( cnst.rel_viol() > feps )
      return( false );
  }
 else
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for(  ; i < get_NStaticArcs() ; ++i ) {
   c_FNumber Ui = get_U( i );
   c_FNumber xi = x[ i ].get_value();
   if( Ui >= Inf<FNumber>() ) {
    if( xi < - feps )
     return( false );
    }
   else {
    c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                               std::max( - xi , xi - Ui ) / std::abs( Ui );
    if( slck > feps )
     return( false );
    }
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for(  ; i < get_NArcs() ; ++i ) {
    c_FNumber Ui = get_U( i );
    c_FNumber xi = (*(dxi++)).get_value();
    if( Ui >= Inf<FNumber>() ) {
     if( xi < - feps )
      return( false );
     }
    else {
     c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                                std::max( - xi , xi - Ui ) / std::abs( Ui );
     if( slck > feps )
      return( false );
     }
    }
   }

 return( true );

 }  // end( MCFBlock::bound_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::dual_feasible( c_CNumber ceps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( E.empty() && dE.empty() )
   throw( std::logic_error( "Constraint required for dual_feasible( , true )"
			    ) );

  if( get_objective().empty() )
   throw( std::logic_error( "Objective required for dual_feasible( , true )"
			    ) );

  auto obj = boost::any_cast<FRealObjective *>( & get_objective() );
  assert( obj );
  #ifdef NDEBUG
   auto lfo = static_cast<const LinearFunction *>( (*obj)->get_function() );
  #else
   auto lfo = dynamic_cast<const LinearFunction *>( (*obj)->get_function() );
   assert( lfo );
  #endif

  for( auto & pi : (*lfo).get_v_var() ) {
   auto xi = pi.first;
   auto RCi = pi.second;
   for( Index j = 0 ; j < xi->get_num_active() ; ++j ) {
    ThinVarDepInterface * ci = xi->get_active( j );
    auto rci = dynamic_cast<FRowConstraint *>( ci );
    if( rci ) {
     #ifdef NDEBUG
      auto lfi = static_cast<const LinearFunction *>( rci->get_function() );
     #else
      auto lfi = dynamic_cast<const LinearFunction *>( rci->get_function() );
      assert( lfi );
     #endif
     RCi -= rci->get_dual() * lfi->get_coefficient( xi );
     }
    else {
     auto bci = dynamic_cast<BoxConstraint *>( ci );
     assert( bci );
     RCi -= bci->get_dual();
     }
    }

   RCi = std::abs( RCi );
   if( pi.second != 0 )
    RCi /= pi.second;

   if( RCi > ceps )
    return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Vec_CNumber RC;
  get_RC( RC );
  Vec_CNumber Pi;
  get_pi( Pi );

  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   c_CNumber Ci = get_C( i );
   c_CNumber RCi = Ci + Pi( SN[ i ] - 1 ) - Pi( EN[ i ] - 1 );
   c_CNumber df = std::abs( RCi - RC[ i ] );
   c_CNumber mx = std::max( std::abs( Ci ) , df );
   if( mx == 0 ) {
    if( df > ceps )
     return( false );
    }
   else
    if( df > ceps * mx )
     return( false );
   }
  }

 return( true );

 }  // end( MCFBlock::dual_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::complementary_slackness( c_CNumber ceps , c_FNumber feps ,
					bool useabstract )
{
 Vec_CNumber RC;
 get_RC( RC );
  
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  assert( ! get_objective().empty() );  // ... which must exist
  auto obj = boost::any_cast<FRealObjective *>( & get_objective() );
  assert( obj );
  #ifdef NDEBUG
  auto lfo = static_cast<const LinearFunction *>( (*obj)->get_function() );
  #else
  auto lfo = dynamic_cast<const LinearFunction *>( (*obj)->get_function() );
   assert( lfo );
  #endif
  Index i = 0;

  // static part
  if( HasStaticX() )
   if( UB.empty() ) {
    for( ; i < get_NStaticArcs() ; ++i ) {
     c_CNumber Ci = lfo->get_coefficient( &x[ i ] );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     if( ( x[ i ]->get_value() > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    }
  else
   for( ; i < get_NStaticArcs() ; ++i ) {
    c_CNumber Ci = lfo->get_coefficient( &x[ i ] );
    CNumber RCi = RC[ i ];
    if( Ci )
     RCi /= Ci;
    c_FNumber xiv = x[ i ]->get_value();
    c_FNumber UBi = UB[ i ].get_rhs();
    if( UBi >= Inf<RowConstraint::RHSValue>() ) {
     if( ( xiv > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    else {
     c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
     if( ( ( xiv > sfeps ) && ( RCi < - ceps ) ) ||
	 ( ( UBi - xiv > sfeps ) && ( RCi > ceps ) ) )
      return( false );
     }
    }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();

   if( dUB.empty() ) {
    for( ; i < get_NStaticArcs() ; ++i ) {
     c_CNumber Ci = lfo->get_coefficient( *dxi );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     if( ( (*(dxi++))->get_value() > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    }
   else {
    auto dubi = dU.begin();

    for( ; i < get_NArcs() ; ++i ) {
     c_CNumber Ci = lfo->get_coefficient( *dxi );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     c_FNumber dxiv = (*(dxi++))->get_value();
     c_FNumber UBi = (*(dubi++)).get_rhs();
     if( UBi >= Inf<RowConstraint::RHSValue>() ) {
      if( ( dxiv > feps ) && ( RCi < - ceps ) )
       return( false );
      }
     else {
      c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
      if( ( ( dxiv > sfeps ) && ( RCi < - ceps ) ) ||
	  ( ( UBi - dxiv > sfeps ) && ( RCi > ceps ) ) )
       return( false );
      }
     }
    }
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for(  ; i < get_NStaticArcs() ; ++i ) {
   c_CNumber Ci = get_C( i );
   CNumber RCi = RC[ i ];
   if( Ci != 0 )
    RCi /= C[ i ];

   c_FNumber Ui = get_U( i );
   c_FNumber xiv = x[ i ].get_value();

   if( Ui >= Inf<FNumber>() ) {
    if( ( xiv > feps ) && ( RCi < - ceps ) )
     return( false );
    }
   else {
    c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
    if( ( ( xiv > sfeps ) && ( RCi < - ceps ) ) ||
	( ( Ui - xiv > sfeps ) && ( RCi > ceps ) ) )
     return( false );
    }
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();

   for(  ; i < get_NArcs() ; ++i ) {
    c_FNumber dxiv = (*(dxi++)).get_value();
    c_CNumber Ci = get_C( i );
    CNumber RCi = RC[ i ];
    if( Ci != 0 )
     RCi /= C[ i ];

    c_FNumber Ui = get_U( i );
    if( Ui >= Inf<FNumber>() ) {
     if( ( dxiv > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    else {
     c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
     if( ( ( dxiv > sfeps ) && ( RCi < - ceps ) ) ||
	 ( ( Ui - dxiv > sfeps ) && ( RCi > ceps ) ) )
      return( false );
     }
    }
   }
  }

 return( true );

 }  // end( MCFBlock::complementary_slackness )

/*--------------------------------------------------------------------------*/

bool MCFBlock::is_feasible( bool useabstract , Configuration *fsbc )
{
 FNumber eps = 0;
 auto tfsbc = dynamic_cast<SimpleConfiguration<FNumber> *>( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast<SimpleConfiguration<FNumber> *>(
			         f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return( flow_feasible( eps , useabstract ) &&
	 bound_feasible( eps , useabstract ) );

 }  //  end( MCFBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::is_optimal( bool useabstract , Configuration *optc )
{
 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  auto toptc =
    dynamic_cast<SimpleConfiguration<std::pair<CNumber,FNumber> > *>( optc );

  if( toptc ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
   }
  else {
   auto ttoptc = dynamic_cast<SimpleConfiguration<CNumber> *>( optc );

   if( ( ! ttoptc ) && f_BlockConfig &&
       f_BlockConfig->f_is_optimal_Configuration )
    ttoptc = dynamic_cast<SimpleConfiguration<CNumber> *>(
			          f_BlockConfig->f_is_optimal_Configuration );
   if( ttoptc )
    ceps = ttoptc->f_value;

   if( f_BlockConfig && f_BlockConfig->f_is_feasible_Configuration ) {
    auto fsbc = dynamic_cast<SimpleConfiguration<FNumber> *>(
			         f_BlockConfig->f_is_feasible_Configuration );
    if( fsbc )
     feps = fsbc->f_value;
    }
   }
  }
 else
  if( f_BlockConfig ) {
   if( f_BlockConfig->f_is_optimal_Configuration ) {
    auto csbc = dynamic_cast<SimpleConfiguration<CNumber> *>(
		                  f_BlockConfig->f_is_optimal_Configuration );
    if( csbc )
     ceps = csbc->f_value;
    }

   if( f_BlockConfig->f_is_feasible_Configuration ) {
    auto fsbc = dynamic_cast<SimpleConfiguration<FNumber> *>(
			        f_BlockConfig->f_is_feasible_Configuration );
    if( fsbc )
     feps = fsbc->f_value;
    }
   }

 return( flow_feasible( feps , useabstract ) &&
	 bound_feasible( feps , useabstract ) &&
	 dual_feasible( ceps , useabstract ) &&
	 complementary_slackness( ceps , feps , useabstract ) );

 }  //  end( MCFBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/

Block * MCFBlock::get_R3_Block( Configuration *r3bc )
{
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 auto MCFB = new MCFBlock();

 MCFB->load( get_NNodes() , EN , SN , U , C , B ,
	     get_NNodes() - get_NStaticNodes() ,
	     get_NArcs() - get_NStaticArcs() ,
	     get_MaxNNodes() - get_NStaticNodes() ,
	     get_MaxNArcs() - get_NStaticArcs() );
 
 return( MCFB );

 }  // end( MCFBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void MCFBlock::map_back_solution( Block *R3B , Configuration *r3bc ,
				               Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast<MCFBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
   tsolc = dynamic_cast<SimpleConfiguration<int> *>(
			            f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map back primal solution - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wsol != 2 ) {
  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( HasStaticX() )
   for( auto xi = x.begin() , r3bxi = MCFB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! xi->is_fixed() )
     xi->set_value( r3bxi->get_value() );
 
  // dynamic part
  // note that if MCFB->dx is longer than this->dx the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for( auto r3bdxi = MCFB->dx.begin() ;
	( dxi != dx.end() ) && ( r3bdxi != MCFB->dx.end() ) ;
	++dxi , ++r3bdxi )
    if( ! dxi->is_fixed() )
     dxi->set_value( r3bdxi->get_value() );

   for( ; dxi != dx.end() ; ++dxi )
    if( ! dxi->is_fixed() )
     dxi->set_value();
   }
  }

 // if required, map back dual solution - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 2 ) && ( E.size() || dE.size() ) ) {  // ... if there is any

  // map back the potentials- - - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticNodes() != get_NStaticNodes() )
   throw( std::invalid_argument( "incompatible static potential size" ) );

  // static part
  if( HasStaticE() )
   for( auto ei = E.begin() , r3bei = MCFB->E.begin() ; ei != E.end() ; )
    (ei++)->set_dual( (r3bei++)->get_dual() );
 
  // dynamic part
  // note that if MCFB->dE is longer than this->dE the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( HasDynamicE() ) {
   auto dei = de.begin();

   for( auto r3bdei = MCFB->dE.begin() ;
	( dei != dE.end() ) && ( r3bdei != MCFB->dE.end() ) ; )
    (dei++)->set_dual( (r3bdei++)->get_dual() );
 
   while( dei != dE.end() )
    (dei++)->set_dual( 0 );
   }

  // map back the reduced costs - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static reduced cost size" ) );

  // static part
  if( HasStaticX() && ( ! UB.empty() ) )
   if( MCFB->UB.empty() ) {
    for( auto & cnst : UB )
     cnst.set_dual( 0 );
    }
   else
    for( auto drci = UB.begin() , r3bdrci = MCFB->UB.begin() ;
	 drci != UB.end() ; )
      (drci++)->set_dual( (r3bdrci++)->get_dual() );

  // dynamic part
  if( HasDynamicX() && ( ! dUB.empty() ) )
   if( MCFB->dUB.empty() ) {
    for( auto & cnst : dUB )
     cnst.set_dual( 0 );
    }
   else {
    auto drci = dUC.begin();
    for( auto r3bdrci = MCFB->dUB.begin() ;
	 ( drci != dUB.end() ) && ( r3bdrci != MCFB->dUB.end() ) ;
	 ++drci , ++r3bdrci )
     drci->set_dual( r3bdrci->get_dual() );
 
    while( drci != dUB.end() )
     (drci++)->set_dual( 0 );
    }
  }
 }  // end( MCFBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void MCFBlock::map_forward_solution( Block *R3B , Configuration *r3bc ,
				                  Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast<MCFBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
  tsolc = dynamic_cast<SimpleConfiguration<int> *>(
			            f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map forward primal solution- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wsol != 2 ) {
  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( MCFB->HasStaticX() )
   for( auto xi = x.begin() , r3bxi = MCFB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! r3bxi->is_fixed() )
     r3bxi->set_value( xi->get_value() );
 
  // dynamic part
  // note that if this->dx is longer than MCFB->dx the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( MCFB->HasDynamicX() ) {
   auto r3bdxi = MCFB->dx.begin();
   for( auto dxi = dx.begin();
	( dxi != dx.end() ) && ( r3bdxi != MCFB->dx.end() ) ;
	++dxi , ++r3bdxi )
    if( ! r3bdxi->is_fixed() )
     r3bdxi->set_value( dxi->get_value() );

   for( ; r3bdxi != MCFB->dx.end() ; ++r3bdxi )
    if( ! r3bdxi->is_fixed() )
     r3bdxi->set_value();
   }
  }

 // if required, map forward dual solution- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 2 ) && ( E.size() || dE.size() ) ) {  // ... if there is any

  // map forward the potentials - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticNodes() != get_NStaticNodes() )
   throw( std::invalid_argument( "incompatible static potential size" ) );

  // static part
  if( MCFB->HasStaticE() )
   for( auto ei = E.begin() , r3bei = MCFB->E.begin() ; ei != E.end() ; )
    (r3bei++)->set_dual( (ei++)->get_dual() );
 
  // dynamic part
  // note that if this->dE is longer than MCFB->dE the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( MCFB->HasDynamicE() ) {
   auto r3bdei = MCFB->dE.begin();

   for( auto dei = de.begin() ;
	( dei != dE.end() ) && ( r3bdei != MCFB->dE.end() ) ; )
    (r3bdei++)->set_dual( (dei++)->get_dual() );
 
   while( r3bdei != MCFB->dE.end() )
    (r3bdei++)->set_dual( 0 );
   }

  // map forward the reduced costs- - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static reduced cost size" ) );

  // static part
  if( MCFB->HasStaticX() && ( ! MCFB->UB.empty() ) )
   if( UB.empty() ) {
    for( auto & cnst : MCFB->UB )
     cnst.set_dual( 0 );
    }
   else
    for( auto drci = UB.begin() , r3bdrci = MCFB->UB.begin() ;
	 drci != UB.end() ; )
      (r3bdrci++)->set_dual( (drci++)->get_dual() );

  // dynamic part
  if( MCFB->HasDynamicX() && ( ! MCFB->dUB.empty() ) )
   if( UB.empty() ) {
    for( auto & cnst : MCFB->dUB )
     cnst.set_dual( 0 );
    }
   else {
    auto r3bdrci = MCFB->dUB.begin();

    for( auto drci = dUB.begin() ;
	 ( drci != dUB.end() ) && ( r3bdrci != MCFB->dUB.end() ) ; )
     (r3bdrci++)->set_dual( (drci++)->get_dual() );
 
    while( r3bdrci != MCFB->dUB.end() )
     (r3bdrci++)->set_dual( 0 );
    }
  }
 }  // end( MCFBlock::map_forward_solution )

/*--------------------------------------------------------------------------*/

bool MCFBlock::map_forward_Modification( Block *R3B , sp_Mod mod ,
					 Configuration *r3bc ,
					 c_ModParam issuePMod ,
					 c_ModParam issueAMod )
{
 auto MCFB = dynamic_cast<MCFBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 /* When a GroupModification is processed, if no channel is provided, then
    one is opened. This only happens "at root", after which in guts_of_mfM()
    whenever a GroupModification is processed, then the channel is nested.
    Indeed, if the "root" Modification is not a GroupModification, then there
    cannot be any GroupModification in it. */

 ModParam iPM = issuePMod;
 ModParam iPA = make_par( std::min( ModParam( eNoBlck ) ,
				    par2mod( issueAMod ) ) ,
			  par2chnl( issueAMod ) );

 const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
 if( tmod ) {  // if the channels are the default ones, open new ones
  if( ! par2chnl( issuePMod ) )
   iPM = make_par( par2concern( issuePMod ) , MCFB->open_channel() );
  if( ! par2chnl( issueAMod ) )
   iPA = make_par( par2concern( issueAMod ) , MCFB->open_channel() );
  }

 /* Use a Lambda to define a "guts" of the method that can be called
    recursively without having to pass "local globals". Note the trick of
    defining the std::function object and "passing" it to the lambda,
    which allows recursive calls. Note the need to explicitly capture
    "this" to use fields/methods of the class. */

 std::function< bool( sp_Mod )> guts_of_mfM;
 guts_of_mfM = [ this , & guts_of_mfM , & MCFB , & iPM , & iPA ]( sp_Mod mod
								  ) {
  // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* This requires to patiently sift through the possible Modification types
     to find what this Modification exactly is, and call the appropriate
     method of either MCFB, for a "physical Modification", or of the "abstract
     representation" of MCFB for an "abstract Modification". */

  //!! std::cout << *mod << std::endl;
  
  // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
   if( tmod ) {
    MCFB->nest_channel( par2chnl( iPM ) );  // nest the channel for PM
    MCFB->nest_channel( par2chnl( iPA ) );  // nest the channel for PA

    bool ok = true;
    for( const auto & submod : tmod->v_sub_Modifications )
     if( ! guts_of_mfM( submod ) )
      ok = false;

    MCFB->un_nest_channel( par2chnl( iPM ) );  // un-nest the channel for PM
    MCFB->un_nest_channel( par2chnl( iPA ) );  // un-nest the channel for PA

    return( ok );
    }
   }

  // MCFBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note: in the following we can assume that C, B and U are nonempty. This
     is because they can be empty only if they are so when the object is
     loaded. But if a Modification has been issued they are no longer empty
     (a Modification changin nothing from the "empty" state is not issued). */
  {
   const auto tmod = std::dynamic_pointer_cast<MCFBlockRngdMod>( mod );
   if( tmod ) {
    switch( tmod->f_type ) {
     case( MCFBlockMod::eChgCost ):
      if( tmod->f_stop == tmod->f_strt + 1 )
       MCFB->chg_cost( C[ tmod->f_strt ] , tmod->f_strt , iPM , iPA );
      else
       MCFB->chg_costs( C.begin() + tmod->f_strt , tmod->f_strt ,
			tmod->f_stop , iPM , iPA );
      break;
     case( MCFBlockMod::eChgCaps ):
      if( tmod->f_stop == tmod->f_strt + 1 )
       MCFB->chg_ucap( U[ tmod->f_strt ] , tmod->f_strt , iPM , iPA );
      else
      MCFB->chg_ucaps( U.begin() + tmod->f_strt , tmod->f_strt ,
		       tmod->f_stop , iPM , iPA );
      break;
     case( MCFBlockMod::eChgDfct ):
      if( tmod->f_stop == tmod->f_strt + 1 )
       MCFB->chg_dfct( B[ tmod->f_strt ] , tmod->f_strt , iPM , iPA );
      else
      MCFB->chg_dfcts( B.begin() + tmod->f_strt , tmod->f_strt ,
		       tmod->f_stop , iPM , iPA );
      break;
     case( MCFBlockMod::eOpenArc ):
      if( tmod->f_stop == tmod->f_strt + 1 )
       MCFB->open_arc( tmod->f_strt , iPM , iPA );
      else
       MCFB->open_arcs( tmod->f_strt , tmod->f_stop , iPM , iPA );
      break;
     case( MCFBlockMod::eCloseArc ):
      if( tmod->f_stop == tmod->f_strt + 1 )
       MCFB->close_arc( tmod->f_strt , iPM , iPA );
      else
      MCFB->close_arcs( tmod->f_strt , tmod->f_stop , iPM , iPA );
      break;
     default:
      throw( std::invalid_argument( "unknown MCFBlockRngdMod type" ) );
     }
    return( true );
    }
   }

  // MCFBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   const auto tmod = std::dynamic_pointer_cast<MCFBlockSbstMod>( mod );
   /* Note that if MCFB will issue a physical modification, then tmod->f_nms
    * need be copied, since then the chg_*() methods "consume" the names
    * vector, but the original vector in mod has to be preserved. Otherwise
    * the names vector is not "consumed", so no copy is needed. */

   if( tmod ) {
    switch( tmod->f_type ) {
     case( MCFBlockMod::eChgCost ): {
      Vec_CNumber NCost( tmod->f_nms.size() );
      for( Index i = 0 ; i < NCost.size() ; i++ )
       NCost[ i ] = C[ tmod->f_nms[ i ] ];

      if( MCFB->issue_pmod( iPM ) )
       MCFB->chg_costs( NCost.begin() ,
			std::move( Vec_Index( tmod->f_nms ) ) , iPM , iPA );
      else
       MCFB->chg_costs( NCost.begin() ,
			std::move( tmod->f_nms ) , iPM , iPA );

      break;
      }
     case( MCFBlockMod::eChgCaps ): {
      Vec_FNumber NCap( tmod->f_nms.size() );
      for( Index i = 0 ; i < NCap.size() ; i++ )
       NCap[ i ] = U[ tmod->f_nms[ i ] ];

      if( MCFB->issue_pmod( iPM ) )
       MCFB->chg_ucaps( NCap.begin() ,
			std::move( Vec_Index( tmod->f_nms ) ) , iPM , iPA );
      else
       MCFB->chg_ucaps( NCap.begin() ,
			std::move( tmod->f_nms ) , iPM , iPA );
      break;
      }
     case( MCFBlockMod::eChgDfct ): {
      Vec_FNumber NDfct( tmod->f_nms.size() );
      for( Index i = 0 ; i < NDfct.size() ; i++ )
       NDfct[ i ] = B[ tmod->f_nms[ i ] ];

      if( MCFB->issue_pmod( iPM ) )
       MCFB->chg_dfcts( NDfct.begin() ,
			std::move( Vec_Index( tmod->f_nms ) ) , iPM , iPA );
      else
       MCFB->chg_dfcts( NDfct.begin() ,
			std::move( tmod->f_nms ) , iPM , iPA );
      break;
      }
     case( MCFBlockMod::eOpenArc ):
      if( MCFB->issue_pmod( iPM ) )
       MCFB->open_arcs( std::move( Vec_Index( tmod->f_nms ) ) , iPM , iPA );
      else
       MCFB->open_arcs( std::move( tmod->f_nms ) , iPM , iPA );
      break;
     case( MCFBlockMod::eCloseArc ):
      if( MCFB->issue_pmod( iPM ) )
       MCFB->close_arcs( std::move( Vec_Index( tmod->f_nms ) ) , iPM , iPA );
      else
       MCFB->close_arcs( std::move( tmod->f_nms ) , iPM , iPA );
      break;
     default:
      throw( std::invalid_argument( "unknown MCFBlockSbstMod type" ) );
     }

    return( true );
    }
   }

  // NBModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   const auto tmod = std::dynamic_pointer_cast<NBModification>( mod );
   if( tmod ) {
    // this is the "nuclear option": the MCFBlock has been re-loaded
    // one should check that the Block is this MCFBlock, but it cannot
    // be otherwise, can it?

    MCFB->load( get_NNodes() , EN , SN , U , C , B );
    return( true );
    }
   }

  return( false );

  };  // end( guts_of_mfM )- - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // finally, call the "guts of"- - - - - - - - - - - - - - - - - - - - - - - -

 bool ok = guts_of_mfM( mod );  // now the actual call

 if( tmod ) {  // now close the opened channels, if any
  if( iPM != issuePMod ) MCFB->close_channel( par2chnl( iPM ) );
  if( iPA != issueAMod ) MCFB->close_channel( par2chnl( iPA ) );
  }

 return( ok );

 }  // end( MCFBlock::map_forward_Modification )

/*--------------------------------------------------------------------------*/

bool MCFBlock::map_back_Modification( Block *R3B , sp_Mod mod ,
				      Configuration *r3bc ,
				      c_ModParam issuePMod ,
				      c_ModParam issueAMod )
{
 /* Fantastically dirty trick: because the two objects are copies, mapping
    back a Modification to this from R3B is the same as mapping forward a
    Modification from R3B to this. */

 auto MCFB = dynamic_cast<MCFBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a MCFBlock" ) );

 return( MCFB->map_forward_Modification( this , mod , r3bc , issuePMod ,
					 issueAMod ) );

 }  // end( MCFBlock::map_back_Modification )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * MCFBlock::get_Solution( Configuration *solc , bool emptys )
{
 int wsol = 0;
 auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
  tsolc = dynamic_cast<SimpleConfiguration<int> *>(
			            f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 auto *sol = new MCFSolution();

 if( wsol != 2 )
  sol->v_x.resize( get_NArcs() );

 if( wsol != 1 )
  sol->v_pi.resize( get_NNodes() );

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( MCFBlock::get_Solution )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_x( Vec_FNumber & FSol , c_Index strt , c_Index stop )
{
 auto FSi = FSol.begin();
 Index i = strt;
 for( ; i < std::min( stop , get_NStaticArcs() ) ; ++i )
  *(FSi++) = x[ i ].get_value();

 if( HasDynamicX() ) {
  auto dxi = dx.begin();
  for( ; i < std::min( stop , get_NArcs() ) ; ++i )
   *(FSi++) = (*(dxi++)).get_value();
  }
 }  // end( MCFBlock::get_x( interval ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_x( Vec_FNumber & FSol , c_Vec_Index & nms )
{
 auto FSi = FSol.begin();
 auto nmsi = nms.begin();

 if( HasDynamicX() ) {
  while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticArcs() ) )
   *(FSi++) = x[ *(nmsi++) ].get_value();

  if( nmsi == nms.end() )
   return;

  auto i = get_NStaticArcs();
  for( auto dxi = dx.begin() ; nmsi != nms.end() ; ++dxi )
   if( *nmsi == i++ ) {
    *(FSi++) = (*dxi).get_value();
    nmsi++;
    }
  }
 else
  while( nmsi != nms.end() ) 
   *(FSi++) = x[ *(nmsi++) ].get_value();

 }  // end( MCFBlock::get_x( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_pi( Vec_CNumber & PSol , c_Index strt , c_Index stop )
{
 if( E.empty() && dE.empty() )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 auto PSi = PSol.begin();
 Index i = strt;
 for( ; i < std::min( stop , get_NStaticNodes() ) ; ++i )
  *(PSi++) = E[ i ].get_dual();

 if( HasDynamicE() ) {
  auto dei = dE.begin();
  for( ; i < std::min( stop , get_NNodes() ) ; ++i )
   *(PSi++) = (*(dei++)).get_dual();
  }
 }  // end( MCFBlock::get_pi( interval ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_pi( Vec_CNumber & PSol , c_Vec_Index & nms )
{
 if( E.empty() && dE.empty() )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 auto PSi = PSol.begin();
 auto nmsi = nms.begin();

 if( HasDynamicE() ) {
  while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticNodes() ) )
   *(PSi++) = E[ *(nmsi++) ].get_dual();

  if( nmsi == nms.end() )
   return;

  auto i = get_NStaticNodes();
  for( auto dei = dE.begin() ; nmsi != nms.end() ; ++dei )
   if( *nmsi == i++ ) {
    *(PSi++) = (*dei).get_dual();
    nmsi++;
    }
  }
 else
  while( nmsi != nms.end() )
   *(PSi++) = E[ *(nmsi++) ].get_dual();

 }  // end( MCFBlock::get_pi( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_rc( Vec_CNumber & RC , c_Index strt , c_Index stop )
{
 if( E.empty() && dE.empty() )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 auto RCi = RC.begin();

 if( UB.empty() && dUB.empty() )
  for( Index i = strt ; i < std::min( stop , get_NArcs() ) ; ++i )
   *(RCi++) = get_C( i ) + get_pi( SN[ i ] - 1 ) - get_pi( EN[ i ] - 1 );
 else {
  Index i = strt;

  if( HasStaticX() )
   for( ; i < std::min( stop , get_NStaticArcs() ) ; ++i )
    *(RCi++) = UB[ i ].get_dual();

  if( HasDynamicX() ) {
   auto dubi = dUB.begin();
   for( ; i < std::min( stop , get_NArcs() ) ; ++i )
    *(RCi++) = (*(dubi++)).get_dual();
   }
  }
 }  // end( MCFBlock::get_rc( interval ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_rc( Vec_CNumber & RC , c_Vec_Index & nms )
{
 if( E.empty() && dE.empty() )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 auto RCi = RC.begin();
 auto nmsi = nms.begin();

 if( UB.empty() && dUB.empty() )
  for( ; nmsi != nms.end() ; ++nmsi ) 
   *(RCi++) = get_C( *nmsi ) + get_pi( SN[ *nmsi ] - 1 )
                             - get_pi( EN[ *nmsi ] - 1 );
 else
  if( HasDynamicX() ) {
   while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticArcs() ) )
     *(RCi++) = UB[ *(nmsi++) ].get_dual();

   if( nmsi == nms.end() )
    return;

   auto i = get_NStaticArcs();
   for( auto dubi = dUB.begin() ; nmsi != nms.end() ; ++dubi )
    if( *nmsi == i++ ) {
     *(RCi++) = (*dubi).get_dual();
     nmsi++;
     }
   }
  else
   while( nmsi != nms.end() ) 
    *(FSi++) = UB[ *(nmsi++) ].get_dual();

 }  // end( MCFBlock::get_rc( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_x( c_Vec_FNumber_it fstrt , c_Vec_FNumber_it fstop ,
		      c_Index strt )
{
 if( std::distance( rcstop , rcstrt ) + strt > get_NArcs() )
  throw( std::invalid_argument( "too many values provided" ) );

 Index i = strt;

 if( HasStaticX() )
  for( auto xi = x.begin() + strt ;
       ( fstrt != fstop ) && ( i < get_NStaticArcs() ) ; ++i )
   (xi++)->set_value( *(fstrt++) );

 if( ( rcstrt == rcstop ) || ( ! HasDynamicX() ) )
  return;

 auto dxi = dx.begin();
 if( i > get_NStaticArcs() )
  dxi = std::next( dxi , i - get_NStaticArcs() );

 while( ; rcstrt < rcstop ; )
  (dxi++)->set_value( *(fstrt++) );

 }  // end( MCFBlock::set_x( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_pi( c_Vec_CNumber_it pstrt , c_Vec_CNumber_it pstop ,
		       c_Index strt )
{
 if( E.empty() && dE.empty() )  // nowhere to put the value in
  return;                       // cowardly (and silently) return

 if( std::distance( pstop , pstrt ) + strt > get_NArcs() )
  throw( std::invalid_argument( "too many values provided" ) );

 Index i = strt;

 if( HasStaticE() )
  for( auto ei = E.begin() + strt ;
       ( pstrt != pstop ) && ( i < get_NStaticArcs() ) ; ++i )
   (ei++)->set_dual( *(pstrt++) );

 if( ( pstrt == pstop ) || ( ! HasDynamicE() ) )
  return;

 auto dei = dE.begin();
 if( i > get_NStaticArcs() )
  dEi = std::next( dEi , i - get_NStaticArcs() );

 while( ; pstrt < pstop ; )
  (dei++)->set_value( *(fpstrt++) );

 }  // end( MCFBlock::set_x( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_rc( c_Vec_CNumber_it rcstrt , c_Vec_CNumber_it rcstop ,
		       c_Index strt )
{
 if( UB.empty() && dUB.empty() )  // nowhere to put the value in
  return;                         // cowardly (and silently) return

 if( std::distance( rcstop , rcstrt ) + strt > get_NArcs() )
  throw( std::invalid_argument( "too many values provided" ) );

 Index i = strt;

 if( HasStaticX() )
  for( auto ubi = UB.begin() + strt ;
       ( rcstrt != rcstop ) && ( i < get_NStaticArcs() ) ; ++i )
   (ubi++)->set_dual( *(rcstrt++) );

 if( ( rcstrt == rcstop ) || ( ! HasDynamicX() ) )
  return;

 auto dubi = dUB.begin();
 if( i > get_NStaticArcs() )
  dubi = std::next( dubi , i - get_NStaticArcs() );

 while( ; rcstrt < rcstop ; )
  (dubi++)->set_dual( *(rcstrt++) );

 }  // end( MCFBlock::set_rc( range ) )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod );
  }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE MCFBlock ---------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::serialize( netCDF::NcGroup && group ) const
{
 group.putAtt( "type" , "MCFBlock" );

 netCDF::NcDim nn = group.addDim( "NNodes" , get_NNodes() );
 netCDF::NcDim na = group.addDim( "NArcs" , get_NArcs() );

 if( get_NNodes() > get_NStaticNodes() )
  group.addDim( "DynNNodes" , get_NNodes() - get_NStaticNodes() );

 if( get_NArcs() > get_NStaticNArcs() )
  group.addDim( "DynNArcs" , get_NArcs() - get_NStaticNArcs() );

 if( get_MaxNNodes() > get_NStaticNodes() )
  group.addDim( "MaxDynNNodes" , get_MaxNNodes() - get_NStaticNodes() );

 if( get_MaxNArcs() > get_NStaticNArcs() )
  group.addDim( "MaxDynNArcs" , get_MaxNArcs() - get_NStaticNArcs() );
 
 std::vector<size_t> startp = { 0 };
 std::vector<size_t> countpa = { get_NArcs() };
 std::vector<size_t> countpn = { get_NNodes() };

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( startp , countpa ,
							    SN.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( startp , countpa ,
							    EN.data() );
 if( C.size() )
  ( group.addVar( "C" , netCDF::NcDouble() , na ) ).putVar( startp , countpa ,
							    C.data() );
 if( U.size() )
  ( group.addVar( "U" , netCDF::NcDouble() , na ) ).putVar( startp , countpa ,
							    U.data() );
 if( B.size() )
  ( group.addVar( "B" , netCDF::NcDouble() , nn ) ).putVar( startp , countpn ,
							    B.data() );
 }  // end( MCFBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::chg_costs( c_Vec_CNumber_it NCost ,
			  c_Index strt , Index stop ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( stop >= get_NArcs() )
  stop = get_NArcs();

 if( ! C.size() ) {
  if( std::all_of( NCost , NCost + ( stop - strt ) ,
		   []( c_CNumber cst ) { return( cst == 0 ); } ) )
   return;

  C.resize( get_MaxNArcs() , 0 );
  }

 if( not_dry_run( issueAMod ) && ( ! get_objective().empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification(s)

  #ifdef NDEBUG
   auto lfo = static_cast<LinearFunction *>( c.get_function() );
  #else
   auto lfo = dynamic_cast<LinearFunction *>( c.get_function() );
   assert( lfo );
  #endif

  Index cnt = 0;
  Vec_CNumber_it cit = C.begin() + strt;
  for( c_Vec_CNumber_it ncit = NCost ;
       ncit < NCost + ( stop - strt ) ; ++ncit , ++cit )
   if( *cit != *ncit ) {
    *cit = *ncit;
    cnt++;  // meanwhile, count how many real changes happen
    }

  if( ! cnt )  // actually nothing has changed
   return;     // avoid the call, hence issuing the abstract Modification

  if( HasDynamicX() ) {  // there are dynamic arcs
   if( stop <= get_NStaticArcs() ) {
    // but all those in the range are static: hence, the range maps
    // into a range of the coefficient, just have to find the extreme
    auto rstrt = lfo->is_active( &x[ strt ] );
    lfo->modify_coefficients( NCost , rstrt , rstrt + ( stop - strt ) ,
			      issueAMod );
    }
   else {
    // the range mixes static and dynamic: the only way is to use the
    // subset version of modify_coefficients()
    LinearFunction::v_coeff_pair pairs;
    auto pi = pairs.begin();

    if( strt < get_NStaticArcs()  )
     for( auto xi = x.begin() + strt ; xi != x.end() ; ) {
      (*pi).first = *(xi++);
      (*(pi++)).second = *(NCost++);
      }

    for( auto dxi = dx.begin() ;
	 dxi != dx.begin() + ( stop - get_NStaticArcs() ) ; ) {
     (*pi).first = *(dxi++);
     (*(pi++)).second = *(NCost++);
     }

    lfo->modify_coefficients( pairs , false , issueAMod );
    }
   }
  else                 // all arcs are static
   lfo->modify_coefficients( NCost , strt , stop , issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCost , NCost + ( stop - strt ) , C.begin() + strt );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				      MCFBlockMod::eChgCost , strt , stop ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::chg_costs( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_costs( c_Vec_CNumber_it NCost , Vec_Index && nms ,
			  const bool ordered  ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( ! C.size() ) {
  if( std::all_of( NCost , NCost + nms.size() ,
		   []( c_CNumber cst ) { return( cst == 0 ); } ) )
   return;

  C.resize( get_MaxNArcs() , 0 );
  }

 if( not_dry_run( issueAMod ) && ( ! get_objective().empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  #ifdef NDEBUG
   auto lfo = static_cast<LinearFunction *>( c.get_function() );
  #else
   auto lfo = dynamic_cast<LinearFunction *>( c.get_function() );
   assert( lfo );
  #endif

  LinearFunction::v_coeff_pair ccp( nms.size() );
  Index cnt = 0;

  if( HasDynamicX() )
   if( ordered ) {
    // static part
    auto nit = nms.begin();
    for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ;
	 ++NCost , ++nit ) {
     if( *nit >= get_NStaticArcs() )
      throw( std::invalid_argument( "invalid arc name" ) );
     if( C[ *nit ] != *NCost ) {
      C[ *nit ] = *NCost;
      ccp[ cnt++ ] = std::make_pair( & x[ *nit ] , *NCost );
      }
     }

    // dynamic part
    auto dxi = dx.begin();
    for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
     if( *(nit++) == i ) {
      if( C[ i ] != *NCost ) {
       C[ i ] = *NCost;
       ccp[ cnt++ ] = std::make_pair( &(*dxi) , *NCost );
       }
      NCost++;
      }
    }
   else {
    // make a vector of pairs < arc index , new cost >
    typedef std::pair< Index , CNumber > index_pair;
    std::vector< index_pair > pairs( nms.size() );
    for( Index i = 0 ; i < nms.size() ; ++i ) {
     if( nms[ i ] >= get_NArcs() )
      throw( std::invalid_argument( "invalid arc name" ) );
     pairs[ i ] = std::make_pair( nms[ i ] , *(NCost++) );
     }

    // sort the vector for increasing index
    std::sort( pairs.begin() , pairs.end ,
	       []( index_pair i , index_pair j )
	       { return( i.first < j.first ); } );

    // static part
    auto pit = pairs.begin();
    for( ; ( pit != pairs.end() ) && ( (*pit).first < get_NStaticArcs() ) ;
	 ++pit )
     if( C[ (*pit).first ] != (*pit).second ) {
      C[ (*pit).first ] = (*pit).second;
      ccp[ cnt++ ] = std::make_pair( & x[ (*pit).first ] , (*pit).second );
      }

    // dynamic part
    auto dxi = dx.begin();
    for( Index i = get_NStaticArcs() ; pairs != pairs.end() ; ++i , ++dxi )
     if( (*pit).first == i ) {
      if( C[ i ] != (*pit).second ) {
       C[ i ] = (*pit).second;
       ccp[ cnt++ ] = std::make_pair( &(*dxi) , (*pit).second );
       }
      pit++;
      }
    }
  else
   for( auto nit = nms.begin() ; nit != nms.end() ; ++NCost , ++nit ) {
    if( C[ *nit ] != *NCost ) {
     if( *nit >= get_NStaticArcs() )
      throw( std::invalid_argument( "invalid arc name" ) );
     C[ *nit ] = *NCost;
     ccp[ cnt++ ] = std::make_pair( & x[ *nit ] , *NCost );
     }
    }

  if( ! cnt )  // actually nothing has changed
   return;     // nothing to change, hence no Modification

  ccp.resize( cnt );

  // note that ccp is ordered if nms was
  lfo->modify_coefficients( ccp , ordered , issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( C , nms , NCost );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  auto mod = std::make_shared<MCFBlockSbstMod>( this ,
                                  MCFBlockMod::eChgCost , std::move( nms ) );

  Block::add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }  // end( MCFBlock::chg_costs( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_cost( c_CNumber NCost , c_Index arc , 
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( ( ! C.size() ) && NCost )
  C.resize( get_MaxNArcs() , 0 );

 if( C[ arc ] == NCost )
  return;

 if( not_dry_run( issueAMod ) && ( ! get_objective().empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  C[ arc ] = NCost;

  #ifdef NDEBUG
   auto lfo = static_cast<LinearFunction *>( c.get_function() );
  #else
   auto lfo = dynamic_cast<LinearFunction *>( c.get_function() );
   assert( lfo );
  #endif

  lfo->modify_coefficient( i2p_x( arc ) , NCost , issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   C[ arc ] = NCost;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				    MCFBlockMod::eChgCost , arc , arc + 1 ) ,
			   Observer::par2chnl( issueMod ) );
 
 }  // end( MCFBlock::chg_cost )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_ucaps( c_Vec_FNumber_it NCap ,
			  c_Index strt , Index stop ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( stop >= get_NArcs() )
  stop = get_NArcs();

 if( stop <= strt )
  return;

 if( ! U.size() ) {
  if( std::all_of( NCap , NCap + nms.size() ,
		   []( c_FNumber cap ) { return( cap >= Inf<FNumber>() ); } ) )
   return;

  U.resize( get_MaxNArcs() , Inf<FNumber>() );
  }

 c_Index ndiff = countdiff( NCap , NCap + ( stop - strt ) , U.begin() + strt );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( ! ( E.empty() && dE.empty() ) ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( UB.empty() && dUB.empty() )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  Index i = strt;

  // static part
  for( ; i < std::min( stop , get_NStaticArcs() ) ;  ++i , ++NCap )
   if( U[ i ] != *NCap ) {
    U[ i ] = *NCap;
    UB[ i ].set_rhs( *NCap , ampar );
    }

  // dynamic part
  for( auto dubi = dUB.begin() ; i < stop ; ++i , ++NCap , ++dubi )
   if( U[ i ] != *NCap ) {
    U[ i ] = *NCap;
    *(dubi++).set_rhs( *NCap , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCap , NCap + ( stop - strt ) , U.begin() + strt );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				      MCFBlockMod::eChgCaps , strt , stop ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::chg_ucaps( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_ucaps( c_Vec_FNumber_it NCap , Vec_Index && nms ,
			  const bool ordered  ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( ! U.size() ) {
  if( std::all_of( NCap , NCap + nms.size() ,
		   []( c_FNumber cap ) { return( cap >= Inf<FNumber>() ); } ) )
   return;

  U.resize( get_MaxNArcs() , Inf<FNumber>() );
  }

 Index ndiff = countdiff( U , nms , NCap , get_NArcs() );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( ! ( E.empty() && dE.empty() ) ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( UB.empty() && dUB.empty() )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  if( HasDynamicX() )
   if( ordered ) {
    // static part
    auto nit = nms.begin();
    for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ;
	 ++NCap , ++nit ) {
     if( U[ *nit ] != *NCap ) {
      U[ *nit ] = *NCap;
      UB[ *nit ].set_rhs( *NCap , ampar );
      }
     }

    // dynamic part
    auto dubi = dUB.begin();
    for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dubi )
     if( *nit == i ) {
      if( U[ i ] != *NCap ) {
       U[ i ] = *NCap;
       (*dubi).set_rhs( *NCap , ampar );
       }
      nit++;
      NCap++;
      }
    }
   else {
    // make a vector of pairs < arc index , new capacity >
    typedef std::pair< Index , FNumber > index_pair;
    std::vector< index_pair > pairs( nms.size() );
    for( Index i = 0 ; i < nms.size() ; ++i )
     pairs[ i ] = std::make_pair( nms[ i ] , *(NCap++) );

    // sort the vector for increasing index
    std::sort( pairs.begin() , pairs.end ,
	       []( index_pair i , index_pair j )
	       { return( i.first < j.first ); } );

    // static part
    auto pit = pairs.begin();
    for( ; ( pit != pairs.end() ) && ( (*pit).first < get_NStaticArcs() ) ;
	 ++pit )
     if( U[ (*pit).first ] != (*pit).second ) {
      U[ (*pit).first ] = (*pit).second;
      UB[ *nit ].set_rhs( *NCap , ampar );
      }

    // dynamic part
    auto dubi = dUB.begin();
    for( Index i = get_NStaticArcs() ; pairs != pairs.end() ; ++i , ++dubi )
     if( (*pit).first == i ) {
      if( U[ i ] != (*pit).second ) {
       U[ i ] = (*pit).second;
       (*dubi).set_rhs( (*pit).second , ampar );
       }
      pit++;
      }
    }
  else
   for( auto nit = nms.begin() ; nit != nms.end() ; ++NCap , ++nit ) {
    if( U[ *nit ] != *NCap ) {
     U[ *nit ] = *NCap;
     UB[ *nit ].set_rhs( *NCap , ampar );
     }
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( U , nms , NCap );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  auto mod = std::make_shared<MCFBlockSbstMod>( this ,
                                  MCFBlockMod::eChgCaps , std::move( nms ) );

  Block::add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }  // end( MCFBlock::chg_ucaps( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_ucap( c_FNumber NCap , c_Index arc ,
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( ( ! U.size() ) && ( NCap < Inf<FNumber>() ) )
  U.resize( get_MaxNArcs() , Inf<FNumber>() );

 if( U[ arc ] == NCap )
  return;

 if( not_dry_run( issueMod ) )
  U[ arc ] = NCap;  // only change the physical representation - - - - - - -

 if( not_dry_run( issueAMod ) && ( ! ( E.empty() && dE.empty() ) ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( UB.empty() && dUB.empty() )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  if( arc < get_NStaticArcs() )
   UB[ arc ].set_rhs( NCap , issueAMod );
  else
   std::next( dUB.begin() , arc - get_NStaticArcs() )->set_rhs( NCap ,
								issueAMod );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				    MCFBlockMod::eChgCaps , arc , arc + 1 ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::chg_ucap )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_dfcts( c_Vec_CNumber_it NDfct ,
			  c_Index strt , Index stop ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( stop >= get_NNodes() )
  stop = get_NNodes();

 if( ! B.size() ) {
  if( std::all_of( NDfct , NDfct + ( stop - strt ) ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.resize( get_MaxNNodes() , 0 );
  }

 c_Index ndiff = countdiff( NDfct , NDfct + ( stop - strt ) ,
			    B.begin() + strt );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( ! ( E.empty() && dE.empty() ) ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( ; i < std::min( stop , get_NStaticNodes() ) ;  ++i , ++NDfct )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    E[ i ].set_both( *NDfct , ampar );
    }

  // dynamic part
  for( auto dei = dE.begin() ; i < stop ; ++i , ++NDfct , ++dei )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    *(dei++).set_both( *NDfct , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NDfct , NDfct + ( stop - strt ) , B.begin() + strt );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				      MCFBlockMod::eChgDfct , strt , stop ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::chg_dfcts( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Vec_Index && nms ,
			  const bool ordered ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( ! B.size() ) {
  if( std::all_of( NDfct , NDfct + nms.size() ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.resize( get_MaxNNodes() , 0 );
  }

 Index ndiff = countdiff( B , nms , NDfct , get_NNodes() );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( ! ( E.empty() && dE.empty() ) ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  if( HasDynamicE() )
   if( ordered ) {
    // static part
    auto nit = nms.begin();
    for( ; ( nit != nms.end() ) && ( *nit < get_NStaticNodes() ) ;
	 ++NDfct , ++nit ) {
     if( B[ *nit ] != *NDfct ) {
      B[ *nit ] = *NDfct;
      E[ *nit ].set_both( *NDfct , ampar );
      }
     }

    // dynamic part
    auto dei = dE.begin();
    for( Index i = get_NStaticNodes() ; nit != nms.end() ; ++i , ++dubi )
     if( *nit == i ) {
      if( B[ i ] != *NDfct ) {
       B[ i ] = *NDfct;
       (*dei).set_both( *NDfct , ampar );
       }
      nit++;
      NDfct++;
      }
    }
   else {
    // make a vector of pairs < arc index , new capacity >
    typedef std::pair< Index , FNumber > index_pair;
    std::vector< index_pair > pairs( nms.size() );
    for( Index i = 0 ; i < nms.size() ; ++i )
     pairs[ i ] = std::make_pair( nms[ i ] , *(NDfct++) );

    // sort the vector for increasing index
    std::sort( pairs.begin() , pairs.end ,
	       []( index_pair i , index_pair j )
	       { return( i.first < j.first ); } );

    // static part
    auto pit = pairs.begin();
    for( ; ( pit != pairs.end() ) && ( (*pit).first < get_NStaticArcs() ) ;
	 ++pit )
     if( B[ (*pit).first ] != (*pit).second ) {
      B[ (*pit).first ] = (*pit).second;
      E[ *nit ].set_both( *NCap , ampar );
      }

    // dynamic part
    auto dei = dE.begin();
    for( Index i = get_NStaticNodes() ; pairs != pairs.end() ; ++i , ++dubi )
     if( (*pit).first == i ) {
      if( B[ i ] != (*pit).second ) {
       B[ i ] = (*pit).second;
       (*dei).set_both( (*pit).second , ampar );
       }
      pit++;
      }
    }
  else
   for( auto nit = nms.begin() ; nit != nms.end() ; ++NDfct , ++nit ) {
    if( B[ *nit ] != *NDfct ) {
     B[ *nit ] = *NDfct;
     E[ *nit ].set_both( *NDfct , ampar );
     }
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( B , nms , NDfct );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared<MCFBlockSbstMod>( this ,
                                 MCFBlockMod::eChgDfct , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }
 }  // end( MCFBlock::chg_dfcts( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::chg_dfct( c_CNumber NDfct , c_Index nde ,
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nde >= get_NNodes() )
  throw( std::invalid_argument( "invalid node name" ) );

 if( ( ! B.size() ) && NDfct )
  B.resize( get_NNodes() , 0 );

 if( B[ nde ] == NDfct )
  return;

 if( not_dry_run( issueMod ) )
  B[ nde ] = NDfct;  // change the physical representation- - - - - - - - - -

 if( not_dry_run( issueAMod ) && ( ! ( E.empty() && dE.empty() ) ) )
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  if( nde < get_NStaticNodes() )
   E[ nde ].set_both( NDfct , issueAMod );
  else
   std::next( dE.begin() , nde - get_NStaticNodes() )->set_both( NDfct ,
								 issueAMod );
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				    MCFBlockMod::eChgDfct , nde , nde + 1 ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::chg_dfct )

/*--------------------------------------------------------------------------*/

void MCFBlock::close_arcs( c_Index strt , Index stop ,
			   c_ModParam issueMod , c_ModParam issueAMod )
{
 if( stop >= get_NArcs() )
  stop = get_NArcs();

 if( stop <= strt )  // nothing to change
  return;            // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;
  Index i = strt;

  // static part
  for( ; i < std::min( stop , get_NStaticArcs() ) ; ++i )
   if( ! x[ i ].is_fixed() )
    ndiff++;

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < stop ; )
   if( ! (*(dxi++)).is_fixed() )
    ndiff++;

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( i = strt ; std::min( stop , get_NStaticArcs() ) ; ++i )
   if( ! x[ i ].is_fixed() ) {
    x[ i ].set_value( 0 );
    x[ i ].is_fixed( true , ampar );
    }

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < stop ; ++dxi )
   if( ! (*dxi).is_fixed() ) {
    (*dxi).set_value( 0 );
    (*dxi).is_fixed( true , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				     MCFBlockMod::eCloseArc , strt , stop ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::close_arcs( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::close_arcs( Vec_Index && nms , const bool ordered  ,
			   c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ; ++nit )
   if( ! x[ *nit ].is_fixed() )
    ndiff++;

  // dynamic part
  auto dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( ! (*dxi).is_fixed() )
     ndiff++;
    ++nit;
    }

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( nit = nms.begin() ; ( nit != nms.end() ) &&
	                   ( *nit < get_NStaticArcs() ) ; ++nit )
   if( ! x[ *nit ].is_fixed() ) {
    x[ *nit ].set_value( 0 );
    x[ *nit ].is_fixed( true , ampar );
    }

  // dynamic part
  dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( ! (*dxi).is_fixed() ) {
     (*dxi).set_value( 0 );
     (*dxi).is_fixed( true , ampar );
     }
    ++nit;
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockSbstMod>( this ,
                                 MCFBlockMod::eCloseArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::close_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::close_arc( c_Index arc ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( xa.is_fixed() )
   return;

  xa.set_value( 0 );

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa.is_fixed( true , issueAMod );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				   MCFBlockMod::eCloseArc , arc , arc + 1 ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::close_arc )

/*--------------------------------------------------------------------------*/

void MCFBlock::open_arcs( c_Index strt , Index stop ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( stop >= get_NArcs() )
  stop = get_NArcs();

 if( stop <= strt )  // nothing to change
  return;            // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
   Index ndiff = 0;
  Index i = strt;

  // static part
  for( ; i < std::min( stop , get_NStaticArcs() ) ; ++i )
   if( x[ i ].is_fixed() )
    ndiff++;

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < stop ; )
   if( (*(dxi++)).is_fixed() )
    ndiff++;

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( i = strt ; std::min( stop , get_NStaticArcs() ) ; ++i )
   if( x[ i ].is_fixed() )
    x[ i ].is_fixed( false , ampar );

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < stop ; ++dxi )
   if( (*dxi).is_fixed() )
    (*dxi).is_fixed( false , ampar );

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				      MCFBlockMod::eOpenArc , strt , stop ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::open_arcs( range ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::open_arcs( Vec_Index && nms , const bool ordered  ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ; ++nit )
   if( x[ *nit ].is_fixed() )
    ndiff++;

  // dynamic part
  auto dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( (*dxi).is_fixed() )
     ndiff++;
    ++nit;
    }

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( nit = nms.begin() ; ( nit != nms.end() ) &&
	                   ( *nit < get_NStaticArcs() ) ; ++nit )
   if( x[ *nit ].is_fixed() )
    x[ *nit ].is_fixed( false , ampar );

  // dynamic part
  dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( (*dxi).is_fixed() )
     (*dxi).is_fixed( false , ampar );
    ++nit;
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockSbstMod>( this ,
                                  MCFBlockMod::eOpenArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::open_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::open_arc( c_Index arc ,
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( ! xa.is_fixed() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa.is_fixed( false , issueAMod );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<MCFBlockRngdMod>( this ,
				    MCFBlockMod::eOpenArc , arc , arc + 1 ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( MCFBlock::open_arc )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MCFBlock::print( std::ostream &output ) const
{
 if( verbosity_lvl != Block::complete ) {  // non-complete version
  // only basic information 
  output << "MCFBlock with: " << NNodes << " nodes and " << SN.size()
	 << " arcs" << std::endl;

  if( verbosity_lvl == Block::high ) {     // print the graph
   for( Index i = 0 ; i < B.size() ; ++i )
    if( B[ i ] != 0 )
     output << "B[ " << i + 1 << " ] = " << B[ i ] << std::endl;

   if( C.size() )
    if( U.size() )
     for( Index i = 0 ; i < C.size() ; ++i ) {
      output << "( " << SN[ i ] << " , " << EN[ i ] << " ): C = " << C[ i ]
	     << ", U = ";
      print_UB( output , U[ i ] );
      output << std::endl;
      }
    else
     for( Index i = 0 ; i < C.size() ; ++i )
      output << "( " << SN[ i ] << " , " << EN[ i ] << " ): C = " << C[ i ]
	     << std::endl;
   else
    if( U.size() )
     for( Index i = 0 ; i < U.size() ; ++i ) {
      output << "( " << SN[ i ] << " , " << EN[ i ] << " ): U = ";
      print_UB( output , U[ i ] );
      output << std::endl;
      }
    else
     output << "all arcs have 0 cost and +Inf upper bound" << std::endl;
   }
  }
 else  {
  // print header in DIMACS standard format
  output << std::endl << "p min " << NNodes << " "<< SN.size() << std::endl;

  // print node descriptors in DIMACS standard format
  for( Index i = 0 ; i < B.size() ; ++i )
   if( B[ i ] != 0 )
    output << "n\t" << i + 1 << "\t" << - B[ i ] << std::endl;

  // print arc descriptors in DIMACS standard format
  if( C.size() )
   if( U.size() )
    for( Index i = 0 ; i < C.size() ; ++i ) {
     output << "a\t" << SN[ i ] << "\t" << EN[ i ] << "\t0\t";
     print_UB( output , U[ i ] );
     output << "\t" << C[ i ] << std::endl;
     }
   else
    for( Index i = 0 ; i < C.size() ; ++i )
     output << "a\t" << SN[ i ] << "\t" << EN[ i ] << "\t0\t+Inf\t"
	    << C[ i ] << std::endl;
  else
   if( U.size() )
    for( Index i = 0 ; i < U.size() ; ++i ) {
     output << "a\t" << SN[ i ] + 1 << "\t" << EN[ i ] + 1 << "\t0\t";
     print_UB( output , U[ i ] );
     output << "\t0" << std::endl;
     }
   else
    for( Index i = 0 ; i < SN.size() ; ++i )
     output << "a\t" << SN[ i ] + 1 << "\t" << EN[ i ] + 1 << "\t0\t+Inf\t0"
	    << std::endl;
  }
 }  // end( MCFBlock::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline MCFBlock::Index MCFBlock::p2i_x( Variable * const var )
{
 auto i = std::distance( x.data() ,
			 static_cast< ColVariable * const >( var ) );
 if( ( i >= 0 ) && ( i < get_NStaticArcs() ) )
  return( i );

 i = get_NStaticArcs();
 for( auto dxi = dx.begin() ; dxi != dx.end() ; ++i , ++dxi )
  if( &(*dxi) == static_cast< ColVariable * >( var ) )
   return( i );

 throw( std::invalid_argument( "invalid arc name" ) );
 return( 0 );
 }

/*--------------------------------------------------------------------------*/

inline ColVariable * MCFBlock::i2p_x( c_Index i )
{
 if( i >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( i < get_NStaticArcs() )
  return( &x[ i ] );
 else
  return( &( *std::next( dx.begin() , i - get_NStaticArcs() ) ) );
 }

/*--------------------------------------------------------------------------*/

inline MCFBlock::Index MCFBlock::p2i_ub( Constraint * const cns )
{
 auto i = std::distance( UB.data() ,
			 static_cast< LB0Constraint * const >( cns ) );
 if( ( i >= 0 ) && ( i < get_NStaticArcs() ) )
  return( i );

 i = get_NStaticArcs();
 for( auto dubi = dUB.begin() ; dubi != dUB.end() ; ++i , ++dubi )
  if( &(*dubi) == static_cast< LB0Constraint * >( cns ) )
   return( i );

 throw( std::invalid_argument( "invalid arc name" ) );
 return( 0 );
 }

/*--------------------------------------------------------------------------*/

inline LB0Constraint * MCFBlock::i2p_ub( c_Index i )
{
 if( i >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( i < get_NStaticArcs() )
  return( &UB[ i ] );
 else
  return( &( *std::next( dUB.begin() , i - get_NStaticArcs() ) ) );
 }

/*--------------------------------------------------------------------------*/

inline MCFBlock::Index MCFBlock::p2i_e( Constraint * const cns )
{
 auto i = std::distance( E.data() ,
			 static_cast< LB0Constraint * const >( cns ) );
 if( ( i >= 0 ) && ( i < get_NStaticNodes() ) )
  return( i );

 i = get_NStaticNodes();
 for( auto dei = dE.begin() ; dei != dE.end() ; ++i , ++dei )
  if( &(*dei) == static_cast< FRowConstraint * >( cns ) )
   return( i );

 throw( std::invalid_argument( "invalid node name" ) );
 return( 0 );
 }

/*--------------------------------------------------------------------------*/

inline FRowConstraint * MCFBlock::i2p_e( c_Index i )
{
 if( i >= get_NNodes() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( i < get_NStaticNodes() )
  return( &E[ i ] );
 else
  return( &( *std::next( dE.begin() , i - get_NStaticArcs() ) ) );
 }

/*--------------------------------------------------------------------------*/

void MCFBlock::guts_of_destructor( void )
{
 /* clear() all Constraint to ensure that they do not bother to un-register
    themselves from Variable that are going to be deleted anyway. Then
    deletes all the "abstract representation", if any. */

 // clear the bound constraints
 for( auto & cnst : UB )
  cnst.clear();
 for( auto & cnst : dUB )
  cnst.clear();

 // clear the flow conservation constraints
 for( auto & cnst : E )
  cnst.clear();
 for( auto & cnst : dE )
  cnst.clear();

 // then delete them all
 dUB.clear();
 UB.clear();
 dE.clear();
 E.clear();

 // clear the objective function
 c.clear();

 // delete all Variable
 dx.clear();
 x.clear();

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer 
 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 }  // end( MCFBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void MCFBlock::guts_of_add_Modification( sp_Mod mod )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
    to find what this Modification exactly is and appropriately mirror the
    changes to the "abstract representation" to the "physical one".

    Note that here we extensively exploit the fact that

        std::distance(  & x[ 0 ] , & x[ i ] ) = i

    because all Variable belong to the same std::vector (array), and similarly
    for the Constraint, in order to efficiently retrieve the index "i" in the
    "phisical representarion" out of pointers in the "abstract
    representation". */

 // TODO: for GroupModification examine the thing in details to recognise
 //       structures, like a bunch of VariableMod corresponding to a set of
 //       arc opening/closures, and react in an optimized way, like with a
 //       single call to [open/close]_arcs() as opposed a single one

 // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
  if( tmod ) {
   for( const auto & submod : tmod->v_sub_Modifications )
    guts_of_add_Modification( submod );
   return;
   }
  }

 // C05FunctionModLin - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  const auto tmod = std::dynamic_pointer_cast<C05FunctionModLin>( mod );
  if( tmod ) {
   if( get_objective().empty() )
    throw( std::invalid_argument( "Modification to non-constructed Objective"
				  ) );

   auto lfo = static_cast<LinearFunction * const>( tmod->f_function );
   if( static_cast<LinearFunction * const>( c.get_function() ) != lfo )
    throw( std::invalid_argument( "Modification to non-Objective" ) );

   if( tmod->v_vars.size() == 1 ) {  // changing just one cost
    const auto var = tmod->v_vars.front();
    const auto cp = lfo->get_v_var();
    c_Index i = p2i_x( var );
    chg_cost( cp[ i ].second , i , eNoBlck , eDryRun );
    }
   else {                            // changing many costs at once
    Vec_CNumber nC( tmod->v_vars.size() );
    Vec_Index nI( tmod->v_vars.size() );
    Index i = 0;

    const auto cp = lfo->get_v_var();

    if( HasDynamicX() ) {
     // use nI as a temporary to store the index of the Variable
     // in the LinearFunction, but it will be overwrittem
     lfo->map_active( tmod->v_vars , nI , true );

     auto dxi = dx.begin();

     if( x.empty() ) {  // there are only dynamic Variable
      for( Index h = 0 ; i < tmod->v_vars.size() ; ++h , ++dxi )
       if( tmod->v_vars[ i ] == &(*dxi) ) {
	nC[ i ] = cp[ nI[ i ] ].second;
	nI[ i++ ] = h;
        }
      }
     else {
      // there are both static and dynamic Variable, which means that
      // x.front() and x.back() are well-defined
      Index h = 0;

      // first part: variables before x.front() (if any)
      for( ; ( i < tmod->v_vars.size() ) &&
	     ( tmod->v_vars[ i ] < &(x.front()) ) ; ++h , ++dxi )
       if( tmod->v_vars[ i ] == &(*dxi) ) {
	nC[ i ] = cp[ nI[ i ] ].second;
	nI[ i++ ] = h;
        }

      // middle part: variables in x (if any)
      for( ; ( i < tmod->v_vars.size() ) &&
	     ( tmod->v_vars[ i ] < &(x.back()) ) ; ++i ) {
       nC[ i ] = cp[ nI[ i ] ].second;
       nI[ i ] = p2i_i( tmod->v_vars[ i ] );
       }

      // last part: variables after x.back() (if any)
      for( ; i < tmod->v_vars.size() ; ++h , ++dxi )
       if( tmod->v_vars[ i ] == &(*dxi) ) {
	nC[ i ] = cp[ nI[ i ] ].second;
	nI[ i++ ] = h;
        }      
      }
     }
    else  // there are only static Variable
     for( const auto & var : tmod->v_vars ) {
      nI[ i ] = p2i_x( var );
      nC[ i ] = cp[ nI[ i ] ].second;
      i++;
      }

    // check if nI is consecutive, if so use the ranged version
    bool cnsctv = true;
    for( auto itnI = nI.begin() ; ; ) {
     auto nxtitnI = ++itnI;
     if( nxtitnI == nI.end() ) break;
     if( *nxtitnI - *itnI != 1 ) { cnsctv = false; break; }
     itnI = nxtitnI;
     }

    if( cnsctv )
     chg_costs( nC.begin() , nI.front() , nI.back() + 1 , eNoBlck , eDryRun );
    else
     chg_costs( nC.begin() , std::move( nI ) , true , eNoBlck , eDryRun );
    }

   return;
   }
  }

 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  const auto tmod = std::dynamic_pointer_cast<RowConstraintMod>( mod );
  if( tmod ) {
   if( ! ( E.empty() && dE.empty() ) )
    throw( std::invalid_argument(
			     "Modification to non-constructed Constraint" ) );

   if( tmod->f_type == RowConstraintMod::eChgRHS ) {
    auto cp = dynamic_cast<LB0Constraint * const>( tmod->f_constraint );
    if( ! cp )
     throw( std::invalid_argument( "Invalid Modification to Constraint" ) );

    chg_ucap( cp->get_rhs() , p2i_ub( cp ) , eNoBlck , eDryRun );
    return;
    }

   if( tmod->f_type == RowConstraintMod::eChgBTS ) {
    auto cp = static_cast<FRowConstraint * const>( tmod->f_constraint );
    if( ! cp )
     throw( std::invalid_argument( "Invalid Modification to Constraint" ) );

    chg_dfct( cp->get_rhs() , p2i_e( cp ) , eNoBlck , eDryRun );
    return;
    }

   throw( std::invalid_argument( "illegal Modification to Constraint" ) );
   }
  }

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  const auto tmod = std::dynamic_pointer_cast<VariableMod>( mod );
  if( tmod ) {
   auto xi = dynamic_cast<ColVariable * const>( tmod->f_variable );
   if( ! xi )
    throw( std::logic_error( "Modification to wrong type of Variable" ) );
   if( ( xi->get_type() != ColVariable::kNonNegative ) &&
       ( xi->get_type() != ColVariable::kNatural ) )
    throw( std::logic_error( "changing type of flow Variable not allowed" ) );
   
   auto i = p2i_x( xi );
   if( xi->is_fixed() )
    close_arc( i , eNoBlck , eDryRun );
   else
    open_arc( i , eNoBlck , eDryRun );

   return;
   }
  }

 throw( std::invalid_argument( "unsupported Modification to MCFBlock" ) );

 }  // end( MCFBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

inline ModParam MCFBlock::make_amod_param( c_ModParam issueAMod ,
					   c_Index num )
{
 if( issue_mod( issueAMod ) ) {
  ChnlName chnl = par2chnl( issueAMod );
  if( num > 1 ) {            // more than one Modification have to be issued
    if( chnl )               // and a channel is provided
     nest_channel( chnl );   // nest the channel
    else                     // sent to default channel
     chnl = open_channel( nullptr , eNoBlck );  /* open a new channel:
                                                 * this creates a
     * GroupModification, which is flagged as "eNoBlck" because it is
     * the "abstract Modification" corresponding to a "physical Modification"
     * already issued and therefore it must not generate any other
     * "physical Modification" */
   }

  return( make_par( eNoBlck , chnl ) );
  }
 else
  return( eNoMod );
 }

/*--------------------------------------------------------------------------*/

inline void MCFBlock::unmake_amod_param( c_ModParam oldiAM ,
					 c_ModParam newiAM , c_Index num )
{
 if( newiAM == eNoMod )
  return;

 ChnlName chnl = par2chnl( newiAM );
 if( num > 1 ) {               // a channel had been opened/nested
  if( par2chnl( oldiAM ) )     // that's "nested"
   un_nest_channel( chnl );    // un-nest it
  else                         // that's "opened"
   close_channel( chnl );      // close it
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF MCFSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void MCFSolution::read( const Block * const block )
{
 auto MCFB = dynamic_cast<const MCFBlock *>( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a MCFBlock" ) );

 // read flows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( v_x.size() < MCFB->get_NArcs() )
  v_x.resize( MCFB->get_NArcs() );

 auto vxi = v_x.begin();

 // static part
 for( auto & xi : MCFB->x )
  *(vxi++) = xi.get_value();

 // dynamic part
 for( auto & xi : MCFB->dx )
  *(vxi++) = xi.get_value();

 // read potentials - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MCFB->E.empty() && MCFB->dE.empty() )  // no potentials available
  return;

 if( v_pi.size() < MCFB->get_NNodes() )
  v_pi.resize( MCFB->get_NNodes() );

 auto vpii = v_pi.begin();
 
 // static part
 for( auto & ei : MCFB->E )
  *(vpii++) = ei.get_dual();

 // dynamic part
 for( auto & ei : MCFB->dE )
  *(vpii++) = ei.get_dual();

 }  // end( MCFSolution::read )

/*--------------------------------------------------------------------------*/

void MCFSolution::write( Block * const block ) 
{
 auto MCFB = dynamic_cast<MCFBlock *>( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a MCFBlock" ) );

 // write flows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( v_x.size() > 0 ) {
  if( v_x.size() < MCFB->get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  auto vxi = v_x.begin();

  // static part
  for( auto & xi : MCFB->x )
   xi.set_value( *(vxi++) );

  // dynamic part
  for( auto dxi = MCFB->dx.begin() ;
       ( dxi != MCFB->dx.end() ) && ( vxi != v_x.end() ) ; )
   (*(dxi++)).set_value( *(vxi++) );
  }

 // write potentials- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( v_pi.empty() )  // no potentials to write
  return;

 if( MCFB->E.empty() && MCFB->dE.empty() )  // no Constraint to write to
  return;

 if( v_pi.size() < MCFB->get_NStaticNodes() )
  throw( std::invalid_argument( "incompatible potential size" ) );

 auto vpii = v_pi.begin();

 // static part
 for( auto & ei : MCFB->E )
  ei.set_dual( *(vpii++) );

 // dynamic part
 for( auto dei = MCFB->dE.begin() ;
      ( dei != MCFB->dE.end() ) && ( vpii != v_pi.end() ) ; )
  (*(dei++)).set_dual( *(vpii++) );

 // write reduced costs (if any)- - - - - - - - - - - - - - - - - - - - - - -
  
 if( MCFB->UB.empty() && MCFB->dUB.empty() )  // no bounds to write to
  return;

 MCFBlock::Index i = 0;

 // static part
 for( auto ubi = MCFB->UB.begin() ; ubi != MCFB->UB.end() ; ++i )
  *(ubi++).set_dual( MCFB->get_C( i ) + v_pi( MCFB->SN[ i ] - 1 )
		                      - v_pi( MCFB->EN[ i ] - 1 ) );
 // dynamic part
 for( auto dubi = MCFB->dUB.begin() ;
      ( dubi != MCFB->dUB.end() ) && ( i < MCFB->get_NArcs() ) ; ++i )
  *(dubi++).set_dual( MCFB->get_C( i ) + v_pi( MCFB->SN[ i ] - 1 )
		                       - v_pi( MCFB->EN[ i ] - 1 ) );
 }  // end( MCFSolution::write )

/*--------------------------------------------------------------------------*/

MCFSolution * MCFSolution::scale( double factor ) const
{
 auto * sol = MCFSolution::clone( true );

 if( v_x.size() > 0 )
  for( MCFBlock::Index i = 0 ; i < v_x.size() ; ++i )
   sol->v_x[ i ] = v_x[ i ] * factor;

 if( v_pi.size() > 0 )
  for( MCFBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   sol->v_pi[ i ] = v_pi[ i ] * factor;

 return( sol );

 }  // end( MCFSolution::scale )

/*--------------------------------------------------------------------------*/

void MCFSolution::sum( const Solution * solution , double multiplier )
{
 auto MCFS = dynamic_cast<const MCFSolution *>( solution );
 if( ! MCFS )
  throw( std::invalid_argument( "solution is not a MCFSolution" ) );

 if( v_x.size() > 0 ) {
  if( v_x.size() != MCFS->v_x.size() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( MCFBlock::Index i = 0 ; i < v_x.size() ; ++i )
   v_x[ i ] = MCFS->v_x[ i ] * multiplier;
  }

 if( v_pi.size() > 0 ) {
  if( v_pi.size() != MCFS->v_pi.size()  )
   throw( std::invalid_argument( "incompatible potential size" ) );

  for( MCFBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   v_pi[ i ] = MCFS->v_pi[ i ] * multiplier;
  }
 }  // end( MCFSolution::sum )

/*--------------------------------------------------------------------------*/

MCFSolution * MCFSolution::clone( bool empty ) const
{
 auto *sol = new MCFSolution();

 if( empty ) {
  if( v_x.size() > 0 )
   sol->v_x.resize( v_x.size() );

  if( v_pi.size() > 0 )
   sol->v_pi.resize( v_pi.size() );
  }
 else {
  sol->v_x = v_x;
  sol->v_pi = v_pi;
  }

 return( sol );

 }  // end( MCFSolution::clone )

/*--------------------------------------------------------------------------*/
/*----------------------- End File MCFBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

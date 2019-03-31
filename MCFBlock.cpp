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
 if( ~ pC.empty() ) {
  C.resize( get_MaxNArcs() , 0 );
  std::copy( pC.begin() , pC.end() , C.begin() );
  }
 if( ~ pU.empty() ) {
  U.resize( get_MaxNArcs() , Inf<FNumber>() );
  std::copy( pU.begin() , pU.end() , U.begin() );
  }
 if( ~ pB.empty() ) {
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
 C.resize( NArcs );
 U.resize( NArcs );
 B.resize( NNodes );

 NStaticNodes = MaxNNodes = NNodes;
 NStaticArcs = NArcs;

 for( auto & el : B )  // all deficits are 0
  el = 0;              // unless otherwise stated

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

 bool reduce = true;
 for( auto cost : C )
  if( cost ) { reduce = false; break; }

 if( reduce )
  C.clear();

 reduce = true;
 for( auto dfct : B )
  if( dfct ) { reduce = false; break; }

 if( reduce )
  B.clear();

 reduce = true;
 for( auto cap : U )
  if( cap < Inf<FNumber>() ) { reduce = false; break; }

 if( reduce )
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

 size_t na = ( group.getDim( "NArcs" ) ).getSize();
 NNodes = ( group.getDim( "NNodes" ) ).getSize();


 NStaticNodes = MaxNNodes = NNodes;
 NStaticArcs = NArcs;


 
 std::vector<size_t> start = { 0 };
 std::vector<size_t> counta = { na };
 std::vector<size_t> countn = { NNodes };

 netCDF::NcVar cst = group.getVar( "C" );
 if( ! cst.isNull() ) {
  C.resize( na );
  cst.getVar( start , counta , C.data() );
  }

 netCDF::NcVar cap = group.getVar( "U" );
 if( ! cap.isNull() ) {
  U.resize( na );
  cap.getVar( start , counta , U.data() );
  }

 netCDF::NcVar dfc = group.getVar( "B" );
 if( ! dfc.isNull() ) {
  B.resize( NNodes );
  dfc.getVar( start , countn , B.data() );
  }

 netCDF::NcVar sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::logic_error( "Starting Nodes not found" ) );

 SN.resize( na );
 sn.getVar( start , counta , SN.data() );

 netCDF::NcVar en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::logic_error( "Ending Nodes not found" ) );

 EN.resize( na );
 en.getVar( start , counta , EN.data() );

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
 if( x.size() == get_NArcs() )  // the variables are there already
  return;                       // nothing to do

 assert( x.size() == 0 );       // this should only happen once

 x.resize( SN.size() );
 for( auto & var : x ) {
  var.is_positive( true , eNoBlck );
  var.set_Block( this );
  }

 add_static_variable( x );

 }  // end( MCFBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void MCFBlock::generate_abstract_constraints( Configuration *stcc )
{
 if( E.size() == get_NNodes() )  // the constraints are there already
  return;                        // nothing to do

 assert( E.size() == 0 );  // this should only happen once

 // generate the node-arc incidence matrix- - - - - - - - - - - - - - - - - -

 E.resize( get_NNodes() );  // ensure E has n constraints

 // each constraint is an equality, i.e., LHS = RHS = B[ i ]
 if( B.size() )
  for( Index i = 0 ; i < get_NNodes() ; ++i )
   E[ i ].set_both( B[ i ] );
 else
  for( auto & ei : E )
   ei.set_both( 0 );

 // number of nonzeroes in each constraint, i.e., #FS( i ) + #BS( i )
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

 // construct the vector of coefficients
 for( Index i = 0 ; i < get_NArcs() ; ++i ) {
  coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( -1 ) );
  coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( 1 ) );
  }

 for( Index i = 0 ; i < get_NNodes() ; ++i ) {
  E[ i ].set_function( new LinearFunction( std::move( coeffs[ i ] ) ,
					   0 , true ) );
  E[ i ].set_Block( this );  // this is done last ==> no Modification
  }

 add_static_constraint( E );

 // generate the bound constraints- - - - - - - - - - - - - - - - - - - - - -
 // if upper bounds are not there, the LB0Constraint can be skipped entirely
 // if the Configuration agrees

 if( U.size() == 0 ) {
  auto tstcc = dynamic_cast<SimpleConfiguration<int> *>( stcc );

  if( ( ! tstcc ) && f_BlockConfig &&
      f_BlockConfig->f_static_constraints_Configuration )
   tstcc = dynamic_cast<SimpleConfiguration<int> *>(
			 f_BlockConfig->f_static_constraints_Configuration );
  if( tstcc && ( tstcc->f_value != 0 ) )
   return;
  }

 auto LU = new std::vector<LB0Constraint>( SN.size() );
 for( Index i = 0 ; i < SN.size() ; ++i ) {
  (*LU)[ i ].set_variable( & x[ i ] , eNoBlck );
  (*LU)[ i ].set_rhs( U[ i ] , eNoBlck );
  (*LU)[ i ].set_Block( this );  // this is done last ==> no Modification
  }

 add_static_constraint( *LU );

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
  if( C.size() )
   for( Index i = 0 ; i < C.size() ; ++i ) {
    p[ i ].first = &x[ i ];
    p[ i ].second = C[ i ];
    }
  else
   for( Index i = 0 ; i < get_NArcs() ; ++i ) {
    p[ i ].first = &x[ i ];
    p[ i ].second = 0;
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
 c.set_function( new LinearFunction( std::move( p ) , 0 , true ) , eNoMod );
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

  assert( E.size() == get_NNodes() );  // ... which must exist

  for( const auto & cnst : E )
   if( cnst.rel_viol() > feps )
    return( false );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Vec_FNumber tB = B;
  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   c_FNumber xi = x[ i ].get_value();
   tB[ SN[ i ] - 1 ] += xi;
   tB[ EN[ i ] - 1 ] -= xi;
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

  assert( E.size() );  // ... which must exist

  if( get_static_constraints().size() > 1 ) {
   auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					   & get_static_constraints()[ 1 ] );
   assert( lbc );
   for( const auto & cnst : **lbc )
    if( cnst.rel_viol() > feps )
     return( false );
   }
  else
   for( const auto & var : x )
    if( ! var.is_feasible() )
     return( false );  
  }
 else
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < SN.size() ; ++i ) {
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

 return( true );

 }  // end( MCFBlock::bound_feasible )

/*--------------------------------------------------------------------------*/

bool MCFBlock::dual_feasible( c_CNumber ceps , bool useabstract )
{
 assert( E.size() );

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
  auto obj_it = lfo->begin();

  for( Index i = 0 , h = 0 ; i < get_NArcs() ; ++i ) {
   CNumber RCi;
   if( static_cast< const ColVariable * >( &(*obj_it) ) == & x[ i ] ) {
    RCi = lfo->get_coefficient( h++ );
    obj_it++;
    }
   else
    RCi = 0;

   c_CNumber Ci = RCi;

   for( Index j = 0 ; j < x[ i ].get_num_active() ; ++j ) {
    ThinVarDepInterface * ci = x[ i ].get_active( j );
    auto rci = dynamic_cast<FRowConstraint *>( ci );
    if( rci ) {
     #ifdef NDEBUG
      auto lfi = static_cast<const LinearFunction *>( rci->get_function() );
     #else
      auto lfi = dynamic_cast<const LinearFunction *>( rci->get_function() );
      assert( lfi );
     #endif
     ThinVarDepInterface::Index pi = lfi->is_active( &x[ i ] );
     assert( pi < Inf<ThinVarDepInterface::Index>() );
     RCi -= rci->get_dual() * lfi->get_coefficient( pi );
     }
    else {
     auto bci = dynamic_cast<BoxConstraint *>( ci );
     assert( bci );
     RCi -= bci->get_dual();
     }
    }

   RCi = std::abs( RCi );
   if( Ci != 0 )
    RCi /= Ci;

   if( RCi > ceps )
    return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < SN.size() ; ++i ) {
   c_CNumber Ci = get_C( i );
   c_CNumber RCi = Ci + E[ SN[ i ] - 1 ].get_dual()
                      - E[ EN[ i ] - 1 ].get_dual();
   c_CNumber df = std::abs( RCi - get_rc( i ) );
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
 assert( E.size() == get_NNodes() );

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
  auto obj_it = lfo->begin();

  auto nnc = boost::any_cast<std::vector<NNConstraint> *>(
					   & get_static_constraints()[ 1 ] );

  auto lbc = nnc ? nullptr : boost::any_cast<std::vector<LB0Constraint> *>(
					   & get_static_constraints()[ 1 ] );


  for( Index i = 0 , h = 0 ; i < SN.size() ; ++i ) {
   CNumber RCi = get_rc( i );
   if( static_cast< const ColVariable * >( &(*obj_it) ) == & x[ i ] ) {
    RCi /= lfo->get_coefficient( h++ );
    obj_it++;
    }

   c_FNumber xi = x[ i ].get_value();
   c_FNumber UBi = nnc ? (**nnc)[ i ].get_rhs() : (**lbc)[ i ].get_rhs();

   if( UBi >= Inf<RowConstraint::RHSValue>() ) {
    if( ( xi > feps ) && ( RCi < - ceps ) )
     return( false );
    }
   else {
    c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
    if( ( ( xi > sfeps ) && ( RCi < - ceps ) ) ||
	( ( UBi - xi > sfeps ) && ( RCi > ceps ) ) )
     return( false );
    }
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < SN.size() ; ++i ) {
   c_CNumber Ci = get_C( i );
   CNumber RCi = get_rc( i );
   if( Ci != 0 )
    RCi /= C[ i ];

   c_FNumber Ui = get_U( i );
   c_FNumber xi = x[ i ].get_value();

   if( Ui >= Inf<FNumber>() ) {
    if( ( xi > feps ) && ( RCi < - ceps ) )
     return( false );
    }
   else {
    c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
    if( ( ( xi > sfeps ) && ( RCi < - ceps ) ) ||
	( ( Ui - xi > sfeps ) && ( RCi > ceps ) ) )
     return( false );
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

 MCFB->load( get_NNodes() , EN , SN , U , C , B );
 
 return( MCFB );

 }  // end( MCFBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void MCFBlock::map_back_solution( Block *R3B , Configuration *r3bc ,
				               Configuration *solc )
{
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

 if( wsol != 2 ) {  // map back primal solution
  if( MCFB->get_NArcs() != get_NArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( Index i = 0 ; i < get_NArcs() ; ++i )
   if( ! x[ i ].is_fixed() )
    x[ i ].set_value( MCFB->x[ i ].get_value() );
  }

 if( wsol != 2 )    // map back dual solution
  if( E.size() ) {  // ... if there is any
   if( MCFB->get_NArcs() != get_NArcs() )
    throw( std::invalid_argument( "incompatible flow size" ) );
   if( MCFB->get_NNodes() != get_NNodes() )
    throw( std::invalid_argument( "incompatible potential size" ) );

   for( Index i = 0 ; i < get_NNodes() ; ++i )
    E[ i ].set_dual( MCFB->E[ i ].get_dual() );

   if( get_static_constraints().size() > 1 ) {
    auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					    & get_static_constraints()[ 1 ] );
    assert( lbc );
    for( Index i = 0 ; i < get_NArcs() ; ++i )
     (**lbc)[ i ].set_dual( MCFB->get_rc( i ) );
    }
   }

}  // end( MCFBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void MCFBlock::map_forward_solution( Block *R3B , Configuration *r3bc ,
				                  Configuration *solc )
{
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

 if( wsol != 2 ) {  // map forward primal solution
  if( MCFB->get_NArcs() != get_NArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( Index i = 0 ; i < get_NArcs() ; ++i )
   if( ! MCFB->x[ i ].is_fixed() )
    MCFB->x[ i ].set_value( x[ i ].get_value() );
  }

 if( wsol != 2 )          // map forward dual solution
  if( MCFB->E.size() ) {  // ... if there is any
   if( MCFB->get_NArcs() != get_NArcs() )
    throw( std::invalid_argument( "incompatible flow size" ) );
   if( MCFB->get_NNodes() != get_NNodes() )
    throw( std::invalid_argument( "incompatible potential size" ) );

   for( Index i = 0 ; i < get_NNodes() ; ++i )
    MCFB->E[ i ].set_dual( E[ i ].get_dual() );

   if( get_static_constraints().size() > 1 ) {
    auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					   & get_static_constraints()[ 1 ] );
    assert( lbc );
    for( Index i = 0 ; i < get_NArcs() ; ++i )
     MCFB->set_rc( (**lbc)[ i ].get_dual() , i );
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

 if( ! emptys)
  sol->read( this );

 return( sol );

 }  // end( MCFBlock::get_Solution )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_x( Vec_FNumber & FSol , c_Index strt , c_Index stp )
{
 for( Index i = 0 ; i < stp - std::min( strt , get_NArcs() ) ; ++i )
  FSol[ i ] = x[ strt + i ].get_value();

 }  // end( MCFBlock::get_x( interval ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_x( Vec_FNumber & FSol , c_Vec_Index & nms )
{
 for( Index i = 0 ; i < nms.size() ; ++i )
  FSol[ i ] = x[ nms[ i ] ].get_value();

 }  // end( MCFBlock::get_x( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_pi( Vec_CNumber & PSol , c_Index strt , c_Index stp )
{
 if( ! E.size() )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 for( Index i = 0 ; i < stp - std::min( strt , get_NNodes() ) ; ++i )
  PSol[ i ] = E[ strt + i ].get_dual();

 }  // end( MCFBlock::get_pi( interval ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_pi( Vec_CNumber & PSol , c_Vec_Index & nms )
{
 if( ! E.size() )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 for( Index i = 0 ; i < nms.size() ; ++i )
  PSol[ i ] = E[ nms[ i ] ].get_dual();

 }  // end( MCFBlock::get_pi( subset ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_rc( Vec_CNumber & RC , c_Index strt , c_Index stp )
{
 if( ! E.size() )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 if( get_static_constraints().size() > 1 ) {
  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					   & get_static_constraints()[ 1 ] );
  assert( lbc );
  for( Index i = 0 ; i < stp - std::min( strt , get_NArcs() ) ; ++i )
   RC[ i ] = (**lbc)[ strt + i ].get_dual();
  }
 else {
  for( Index i = 0 ; i < stp - std::min( strt , get_NArcs() ) ; ++i )
   RC[ i ] = get_C( strt + i ) + E[ SN[ strt + i ] - 1 ].get_dual()
                               - E[ EN[ strt + i ] - 1 ].get_dual();
  }
 }  // end( MCFBlock::get_rc( interval ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::get_rc( Vec_CNumber & RC , c_Vec_Index & nms )
{
 if( ! E.size() )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 if( get_static_constraints().size() > 1 ) {
  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					   & get_static_constraints()[ 1 ] );
  assert( lbc );
  for( Index i = 0 ; i < nms.size() ; ++i )
   RC[ i ] = (**lbc)[ nms[ i ] ].get_dual();
  }
 else
  for( Index i = 0 ; i < nms.size() ; ++i )
   RC[ i ] = get_C( nms[ i ] ) + E[ SN[ nms[ i ] ] - 1 ].get_dual()
                               - E[ EN[ nms[ i ] ] - 1 ].get_dual();

 }  // end( MCFBlock::get_rc( subset ) )

/*--------------------------------------------------------------------------*/

MCFBlock::CNumber MCFBlock::get_rc( c_Index arc )
{
 if( ! E.size() )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't"
			   ) );
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( get_static_constraints().size() > 1 ) {
  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					     get_static_constraints()[ 1 ] );
  assert( lbc );
  return( (*lbc)[ arc ].get_dual() );
  }
 else
  return( get_C( arc ) + E[ SN[ arc ] - 1 ].get_dual()
	               - E[ EN[ arc ] - 1 ].get_dual() );
 
 }  // end( MCFBlock::get_rc( one ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_rc( c_Vec_CNumber_it rcstrt , c_Vec_CNumber_it rcstop ,
		       c_Index strt )
{
 if( ( ! E.size() ) || ( get_static_constraints().size() <= 1 ) )
  return;  // nowhere to put the value:cowardly (and silently) return

 if( std::distance( rcstop , rcstrt ) + strt > get_NArcs() )
  throw( std::invalid_argument( "too many values provided" ) );

 auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					     get_static_constraints()[ 1 ] );
 assert( lbc );
 for( auto bi = lbc->begin() + strt ; rcstrt < rcstop ; )
  (bi++)->set_dual( *(rcstrt++) );

 }  // end( MCFBlock::set_rc( one ) )

/*--------------------------------------------------------------------------*/

void MCFBlock::set_rc( c_CNumber RC , c_Index arc )
{
 if( ( ! E.size() ) || ( get_static_constraints().size() <= 1 ) )
  return;  // nowhere to put the value:cowardly (and silently) return

 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					     get_static_constraints()[ 1 ] );
 assert( lbc );
 (*lbc)[ arc ].set_dual( RC );

 }  // end( MCFBlock::set_rc( one ) )

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

 std::vector<size_t> startp = { 0 };
 std::vector<size_t> countpa = { get_NArcs() };
 std::vector<size_t> countpn = { get_NNodes() };

 if( C.size() )
  ( group.addVar( "C" , netCDF::NcDouble() , na ) ).putVar( startp , countpa ,
							    C.data() );
 if( U.size() )
  ( group.addVar( "U" , netCDF::NcDouble() , na ) ).putVar( startp , countpa ,
							    U.data() );
 if( B.size() )
  ( group.addVar( "B" , netCDF::NcDouble() , nn ) ).putVar( startp , countpn ,
							    B.data() );

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( startp , countpa ,
							    SN.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( startp , countpa ,
							    EN.data() );

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

 c_Vec_CNumber_it ncit = NCost;
 const c_Vec_CNumber_it ncstp = ncit + ( stop - strt );

 if( ! C.size() ) {
  for( ; ncit < ncstp ; ++ncit )
   if( *ncit )
    break;

  if( ncit >= ncstp )
   return;

  C.resize( get_NArcs() , 0 );
  ncit = NCost;
  }

 Vec_CNumber_it cit = C.begin() + strt;

 if( not_dry_run( issueAMod ) && ( ! get_objective().empty()  ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification(s)

  #ifdef NDEBUG
   auto lfo = static_cast<LinearFunction *>( c.get_function() );
  #else
   auto lfo = dynamic_cast<LinearFunction *>( c.get_function() );
   assert( lfo );
  #endif

  /*!! if( lfo->get_num_active_var() == get_NArcs() ) !!*/ {
   // "dense" objective
   Index cnt = 0;
   for( ; ncit < ncstp ; ++ncit , ++cit )
    if( *cit != *ncit ) {
     *cit = *ncit;
     cnt++;  // meanwhile, count how many real changes happen
     }

   if( ! cnt )  // actually nothing has changed
    return;     // avoid the call, hence issuing the abstract Modification

   lfo->modify_coefficients( NCost , strt , stop , issueAMod );
   }
  /*!!
  else {                                            // "sparse" objective
   Index addv = 0;  // Variable to be added
   Index rmvv = 0;  // Variable to be removed
   Index chgv = 0;  // coefficients to be changed
 
   Vec_CNumber_it cit = C.begin() + strt;
   for( ; ncit < ncstp ; ++ncit , ++cit )
    if( *ncit != *cit ) {
     if( *ncit == 0 )
      rmvv++;
     else
      if( *cit == 0 )
       addv++;
      else
       chgv++;
     }

   if( ! ( addv + rmvv + chgv ) )
    return;

   LinearFunction::v_coeff_pair acp( addv );
   LinearFunction::v_coeff_pair ccp( chgv );
   Vec_p_Var rcp( rmvv );

   addv = 0;
   rmvv = 0;
   chgv = 0;

   // compute the three sets of removed, added and changed coefficients,
   // all the while doing the change, so as to ensure that the change is
   // in place the moment the Modification is issued.
   auto xit = x.begin();
   for( cit = C.begin() + strt , ncit = NCost ; ncit < ncstp ;
	++ncit , ++cit , ++xit )
    if( *ncit != *cit ) {
     ColVariable * xi = & (*xit);
     if( *ncit == 0 )
      rcp[ rmvv++ ] = xi;
     else
      if( *cit == 0 )
       acp[ addv++ ] = std::make_pair( xi , *ncit );
      else
       ccp[ chgv++ ] = std::make_pair( xi , *ncit );

     *cit = *ncit;
     }

   c_Index nmod = ( addv > 0 ) + ( rmvv > 0 ) + ( chgv > 0 );
   c_ModParam ampar = make_amod_param( issueAMod , nmod );

   // note that the vectors are ordered by construction
   if( rmvv ) lfo->remove_variables( std::move( rcp ) , true , ampar );
   if( addv ) lfo->add_variables( std::move( acp ) , true , ampar );
   if( chgv ) lfo->modify_coefficients( std::move( ccp ) , true , ampar );

   unmake_amod_param( issueAMod , ampar , nmod );
   }
   !!*/
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   while( ncit < ncstp )
    *(cit++) = *(ncit++);

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
 assert( nms.size() <= get_NArcs() );

 c_Vec_CNumber_it ncit = NCost;
 const c_Vec_CNumber_it ncstp = ncit + nms.size();

 if( ! C.size() ) {
  for( ; ncit < ncstp ; ++ncit )
   if( *ncit )
    break;

  if( ncit >= ncstp )
   return;

  C.resize( get_NArcs() , 0 );
  ncit = NCost;
  }

 c_Vec_Index_it nit = nms.begin();

 if( not_dry_run( issueAMod ) && ( ! get_objective().empty()  ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  #ifdef NDEBUG
   auto lfo = static_cast<LinearFunction *>( c.get_function() );
  #else
   auto lfo = dynamic_cast<LinearFunction *>( c.get_function() );
   assert( lfo );
  #endif

  /*!! if( lfo->get_num_active_var() == get_NArcs() ) !!*/ {
   // "dense" objective
   LinearFunction::v_coeff_pair ccp( nms.size() );

   Index cnt = 0;
   for( ; ncit < ncstp ; ++ncit , ++nit )
    if( C[ *nit ] != *ncit ) {
     C[ *nit ] = *ncit;
     ccp[ cnt++ ] = std::make_pair( & x[ *nit ] , *ncit );
     }

   if( ! cnt )  // actually nothing has changed
    return;     // nothing to change, hence no Modification

   ccp.resize( cnt );

   // note that ccp is ordered if nms was
   lfo->modify_coefficients( ccp , ordered , issueAMod );
   }
  /*!!
  else {                                            // "sparse" objective
   Index addv = 0;  // Variable to be added
   Index rmvv = 0;  // Variable to be removed
   Index chgv = 0;  // coefficients to be changed

   for( ; ncit < ncstp ; ++ncit , ++nit ) {
    #ifndef NDEBUG
     if( *nit >= get_NArcs() )
      throw( std::invalid_argument( "invalid arc name" ) );
    #endif
    if( *ncit != C[ *nit ] ) {
     if( *ncit == 0 )
      rmvv++;
     else
      if( C[ *nit ] == 0 )
       addv++;
      else
       chgv++;
     }
    }

   if( ! ( addv + rmvv + chgv ) )
    return;

   LinearFunction::v_coeff_pair acp( addv );
   LinearFunction::v_coeff_pair ccp( chgv );
   Vec_p_Var rcp( rmvv );

   addv = 0;
   rmvv = 0;
   chgv = 0;

   // compute the three sets of removed, added and changed coefficients,
   // all the while doing the change, so as to ensure that the change is
   // in place the moment the Modification is issued
   for( nit = nms.begin() , ncit = NCost ; ncit < ncstp ; ++ncit , ++nit )
    if( C[ *nit ] != *ncit ) {
     ColVariable * xi = & x[ *nit ];
     if( *ncit == 0 )
      rcp[ rmvv++ ] = xi;
     else
      if( C[ *nit ] == 0 )
       acp[ addv++ ] = std::make_pair( xi , *ncit );
      else
       ccp[ chgv++ ] = std::make_pair( xi , *ncit );

     C[ *nit ] = *ncit;
     }

   c_Index nmod = ( addv > 0 ) + ( rmvv > 0 ) + ( chgv > 0 );
   c_ModParam ampar = make_amod_param( issueAMod , nmod );

   // note that the vectors are ordered if nms is
   if( rmvv ) lfo->remove_variables( std::move( rcp ) , ordered , nmod );
   if( addv ) lfo->add_variables( std::move( acp ) , ordered , nmod );
   if( chgv ) lfo->modify_coefficients( std::move( ccp ) , ordered , nmod );

   unmake_amod_param( issueAMod , ampar , nmod );
   }
   !!*/
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   for( ; ncit < ncstp ; ++ncit , ++nit )
    C[ *nit ] = *ncit;

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
   C.resize( get_NArcs() , 0 );

 if( C[ arc ] == NCost )
  return;

 if( not_dry_run( issueAMod ) && ( ! get_objective().empty()  ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  #ifdef NDEBUG
   auto lfo = static_cast<LinearFunction *>( c.get_function() );
  #else
   auto lfo = dynamic_cast<LinearFunction *>( c.get_function() );
   assert( lfo );
  #endif

  /*!! if( lfo->get_num_active_var() == get_NArcs() ) !!*/ {
   // "dense" objective
   C[ arc ] = NCost;                                // modify coefficient
   lfo->modify_coefficient( &x[ arc ] , NCost , issueAMod );
   }
  /*!!
  else                                              // "sparse" objective
   if( NCost == 0 ) {        // remove one Variable
    C[ arc ] = 0;
    lfo->remove_variable( &x[ arc ] , issueAMod );
    }
   else
    if( C[ arc ] == 0 ) {    // add one Variable
     C[ arc ] = NCost;
     lfo->add_variable( &x[ arc ] , NCost , issueAMod );
     }
    else {                   // modify one coefficient
     C[ arc ] = NCost;
     lfo->modify_coefficient( &x[ arc ] , NCost , issueAMod );
     }
     !!*/
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

 c_Vec_FNumber_it ncit = NCap;
 const c_Vec_CNumber_it ncstp = ncit + ( stop - strt );

 if( ! U.size() ) {
  for( ; ncit < ncstp ; ++ncit )
   if( *ncit < Inf<FNumber>() )
    break;

  if( ncit >= ncstp )
   return;

  U.resize( get_NArcs() , Inf<FNumber>() );
  ncit = NCap;
  }

 Vec_FNumber_it uit = U.begin() + strt;

 Index ndiff = 0;
 while( ncit < ncstp )
  if( *(ncit++) != *(uit++) )
   ndiff++;

 if( ! ndiff )
  return;

 ncit = NCap;
 uit = U.begin() + strt;

 if( not_dry_run( issueAMod ) && ( ! E.empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( get_static_constraints().size() <= 1 )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					  & get_static_constraints()[ 1 ] );
  assert( lbc );

  for( auto lbit = (*lbc)->begin() + strt ; ncit < ncstp ;
       ++ncit , ++uit , ++lbit )
   if( *uit != *ncit ) {
    *uit = *ncit;
    lbit->set_rhs( *ncit , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   while( ncit < ncstp )
    *(uit++) = *(ncit++);

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
 assert( nms.size() <= get_NArcs() );

 c_Vec_FNumber_it ncit = NCap;
 const c_Vec_FNumber_it ncstp = ncit + nms.size();

 if( ! U.size() ) {
  for( ; ncit < ncstp ; ++ncit )
   if( *ncit < Inf<FNumber>() )
    break;

  if( ncit >= ncstp )
   return;

  U.resize( get_NArcs() , Inf<FNumber>() );
  ncit = NCap;
  }

 c_Vec_Index_it nit = nms.begin();

 Index ndiff = 0;
 for( ; ncit < ncstp ; ++ncit , ++nit ) {
  #ifndef NDEBUG
   if( *nit >= get_NArcs() )
    throw( std::invalid_argument( "invalid arc name" ) );
  #endif
  if( *ncit != U[ *nit ] )
   ndiff++;
  }

 if( ! ndiff )
  return;

 ncit = NCap;
 nit = nms.begin();

 if( not_dry_run( issueAMod ) && ( ! E.empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( get_static_constraints().size() <= 1 )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					  & get_static_constraints()[ 1 ] );
  assert( lbc );

  for( ; ncit < ncstp ; ++ncit , ++nit ) {
   if( U[ *nit ] != *ncit ) {
    U[ *nit ] = *ncit;
    (**lbc)[ *nit ].set_rhs( *ncit , ampar );
    }
   }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   while( ncit < ncstp )
    U[ *(nit++) ] = *(ncit++);

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
  U.resize( get_NArcs() , Inf<FNumber>() );

 if( U[ arc ] == NCap )
  return;

 if( not_dry_run( issueMod ) )
  U[ arc ] = NCap;  // only change the physical representation - - - - - - -

 if( not_dry_run( issueAMod ) && ( ! E.empty() ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( get_static_constraints().size() <= 1 )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					  & get_static_constraints()[ 1 ] );
  assert( lbc );

  (**lbc)[ arc ].set_rhs( NCap , issueAMod );
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

 c_Vec_FNumber_it ndit = NDfct;
 const c_Vec_CNumber_it ndstp = ndit + ( stop - strt );

 if( ! B.size() ) {
  for( ; ndit < ndstp ; ++ndit )
   if( *ndit )
    break;

  if( ndit >= ndstp )
   return;

  B.resize( get_NNodes() , 0 );
  ndit = NDfct;
  }

 Vec_FNumber_it bit = B.begin() + strt;

 Index ndiff = 0;
 while( ndit < ndstp )
  if( *(ndit++) != *(bit++) )
   ndiff++;

 if( ! ndiff )
  return;

 ndit = NDfct;
 bit = B.begin() + strt;

 if( not_dry_run( issueAMod ) && ( ! E.empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  for( auto eit = E.begin() ; ndit < ndstp ; ++ndit , ++bit , ++eit )
   if( *bit != *ndit ) {
    *bit = *ndit;
    eit->set_both( *ndit , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   while( ndit < ndstp )
    *(bit++) = *(ndit++);

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
 assert( nms.size() <= get_NNodes() );

 c_Vec_FNumber_it ndit = NDfct;
 const c_Vec_FNumber_it ndstp = ndit + nms.size();

 if( ! B.size() ) {
  for( ; ndit < ndstp ; ++ndit )
   if( *ndit )
    break;

  if( ndit >= ndstp )
   return;

  B.resize( get_NNodes() , 0 );
  ndit = NDfct;
  }

 c_Vec_Index_it nit = nms.begin();

 Index ndiff = 0;
 for( ; ndit < ndstp ; ++ndit , ++nit ) {
  #ifndef NDEBUG
   if( *nit >= get_NNodes() )
    throw( std::invalid_argument( "invalid node name" ) );
  #endif
  if( *ndit != B[ *nit ] )
   ndiff++;
  }

 if( ! ndiff )
  return;

 ndit = NDfct;
 nit = nms.begin();

 if( not_dry_run( issueAMod ) && ( ! E.empty() ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  for( ; ndit < ndstp ; ++ndit , ++nit )
   if( B[ *nit ] != *ndit ) {
    B[ *nit ] = *ndit;
    E[ *nit ].set_both( *ndit , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   while( ndit < ndstp )
    B[ *(nit++) ] = *(ndit++);

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

 if( not_dry_run( issueAMod ) && ( ! E.empty() ) )
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  E[ nde ].set_both( NDfct , issueAMod );

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
  for( Index i = strt ; i < stop ; ++i )
   if( ! x[ i ].is_fixed() )
    ndiff++;

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  for( Index i = strt ; i < stop ; ++i )
   if( ! x[ i ].is_fixed() ) {
    x[ i ].set_value( 0 );
    x[ i ].is_fixed( true , ampar );
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
 assert( nms.size() <= get_NArcs() );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;
  for( auto i : nms ) {
   #ifndef NDEBUG
    if( i >= get_NArcs() )
     throw( std::invalid_argument( "invalid arc name" ) );
   #endif
   if( ! x[ i ].is_fixed() )
    ndiff++;
   }

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  for( auto i : nms )
   if( ! x[ i ].is_fixed() ) {
    x[ i ].set_value( 0 );
    x[ i ].is_fixed( true , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  auto mod = std::make_shared<MCFBlockSbstMod>( this ,
                                 MCFBlockMod::eCloseArc , std::move( nms ) );

  Block::add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
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
  if( x[ arc ].is_fixed() )
   return;

  x[ arc ].set_value( 0 );

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  x[ arc ].is_fixed( true , issueAMod );
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
  for( Index i = strt ; i < stop ; ++i )
   if( x[ i ].is_fixed() )
    ndiff++;

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  for( Index i = strt ; i < stop ; ++i )
   if( x[ i ].is_fixed() )
    x[ i ].is_fixed( false , ampar );

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
 assert( nms.size() <= get_NArcs() );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;
  for( auto i : nms ) {
   #ifndef NDEBUG
    if( i >= get_NArcs() )
     throw( std::invalid_argument( "invalid arc name" ) );
   #endif
   if( x[ i ].is_fixed() )
    ndiff++;
   }

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  for( auto i : nms )
   if( x[ i ].is_fixed() )
    x[ i ].is_fixed( false , ampar );

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  auto mod = std::make_shared<MCFBlockSbstMod>( this ,
                                   MCFBlockMod::eOpenArc , std::move( nms ) );

  Block::add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
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
  if( ! x[ arc ].is_fixed() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  x[ arc ].is_fixed( false , issueAMod );
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

inline MCFBlock::Index MCFBlock::p2i( Variable * const var )
{
 return( std::distance( x.data() , static_cast< ColVariable * const >( var ) )
	 );
 }

/*--------------------------------------------------------------------------*/

void MCFBlock::guts_of_destructor( void )
{
 /* clear all Constraint to ensure that they do not make any reference to
    no-longer-existing Variable while they are destroyed. */

 assert( ( get_static_constraints().size() == 0 ) ||
	 ( get_static_constraints().size() == 2 ) );

 if( get_static_constraints().size() > 1 ) { // delete LB0 constraints, if any
  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					   & get_static_constraints()[ 1 ] );
  assert( lbc );
  for( auto & cnst : **lbc )  // first clear all the Constraint
   cnst.clear();
  delete *lbc;                // then delete them
  }

 for( Index i = E.size() ; i-- ; )
  E[ i ].clear();

 E.clear();

 // clear the objective function
 c.clear();

 // explicitly clear static Constraint and Variable
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 x.clear();

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
    c_Index i = p2i( var );
    chg_cost( cp[ i ].second , i , eNoBlck , eDryRun );
    }
   else {                            // changing many costs at once
    Vec_CNumber nC( tmod->v_vars.size() );
    Vec_Index nI( tmod->v_vars.size() );
    Index i = 0;
   
    const auto cp = lfo->get_v_var();
    for( const auto & var : tmod->v_vars ) {
     nI[ i ] = p2i( var );
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
   if( ! E.size() )
    throw( std::invalid_argument( "Modification to non-constructed Constraint"
				  ) );

   if( tmod->f_type == RowConstraintMod::eChgRHS ) {
    auto cp = dynamic_cast<LB0Constraint * const>( tmod->f_constraint );
    if( ! cp )
     throw( std::invalid_argument( "Invalid Modification to Constraint" ) );
     
    auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					  & get_static_constraints()[ 1 ] );
    assert( lbc );

    auto i = std::distance( &((*lbc)->front()) , cp );
    if( ( i < 0 ) || ( i >= get_NArcs() ) )
     throw( std::invalid_argument(
			    "Modification to Constraint of another Block" ) );

    chg_ucap( cp->get_rhs() , i , eNoBlck , eDryRun );
    return;
    }

   if( tmod->f_type == RowConstraintMod::eChgBTS ) {
    auto cp = static_cast<FRowConstraint * const>( tmod->f_constraint );
    if( ! cp )
     throw( std::invalid_argument( "Invalid Modification to Constraint" ) );

    auto i = std::distance( &(E.front()) , cp );
    if( ( i < 0 ) || ( i >= NNodes ) )
     throw( std::invalid_argument(
			    "Modification to Constraint of another Block" ) );

    chg_dfct( cp->get_rhs() , i , eNoBlck , eDryRun );
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
   
   auto i = std::distance( &(x.front()) , xi );
   if( ( i < 0 ) || ( i >= get_NArcs() ) )
    throw( std::invalid_argument(
			     "Modification to Variable of another Block" ) );
   if( xi->is_fixed() )
    close_arc( i ,  eNoBlck , eDryRun );
   else
    open_arc( i ,  eNoBlck , eDryRun );

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

 if( v_x.size() > 0 ) {
  if( v_x.size() != MCFB->get_NArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( MCFBlock::Index i = 0 ; i < v_x.size() ; ++i )
   v_x[ i ] = MCFB->x[ i ].get_value();
  }

 if( v_pi.size() > 0 ) {
  if( v_pi.size() != MCFB->get_NNodes() )
   throw( std::invalid_argument( "incompatible potential size" ) );

  if( ! MCFB->E.size() )
   throw( std::invalid_argument( "potential solution not available" ) );

  for( MCFBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   v_pi[ i ] = MCFB->E[ i ].get_dual();
  }
 }  // end( MCFSolution::read )

/*--------------------------------------------------------------------------*/

void MCFSolution::write( Block * const block ) 
{
 auto MCFB = dynamic_cast<MCFBlock *>( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a MCFBlock" ) );

 if( v_x.size() > 0 ) {
  if( v_x.size() != MCFB->get_NArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( MCFBlock::Index i = 0 ; i < v_x.size() ; ++i )
   if( ! MCFB->x[ i ].is_fixed() )
    MCFB->x[ i ].set_value( v_x[ i ] );
  }

 if( v_pi.size() > 0 ) {
  if( v_pi.size() != MCFB->get_NNodes() )
   throw( std::invalid_argument( "incompatible potential size" ) );

  if( ! MCFB->E.size() )
   throw( std::invalid_argument( "potential solution not available" ) );

  for( MCFBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   MCFB->E[ i ].set_dual( v_pi[ i ] );

  if( MCFB->get_static_constraints().size() > 1 ) {
   auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
				     & MCFB->get_static_constraints()[ 1 ] );
   assert( lbc );
   for( MCFBlock::Index i = 0 ; i < MCFB->get_NArcs() ; ++i )
    (**lbc)[ i ].set_dual( MCFB->C[ i ] - v_pi[ MCFB->SN[ i ] ]
		                        + v_pi[ MCFB->EN[ i ] ] );
   }
  }
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

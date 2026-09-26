/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tester of the MCFBlock alone, i.e., with nothing but the core SMS++
 * library: the instances are built in memory, and no Solver of the MCF is
 * involved, the only one attached being a FakeSolver collecting the
 * Modification the MCFBlock issues.
 *
 * The DIMACS and netCDF round trips are checked on instances with no arcs,
 * with unbalanced deficits, with infinite capacities, with closed and with
 * deleted arcs. Every change of the physical data (costs, capacities,
 * deficits, closing and opening arcs, adding and removing arcs), in the
 * single, Range and Subset forms, empty ones included, is checked against
 * the abstract representation, which has to follow it, and against the
 * Modification the FakeSolver receives, under the eNoMod and eModBlck
 * values of the parameters too; so is the other way round, a change of the
 * abstract representation that the MCFBlock has to translate into its
 * data. A flow set by hand in the Variable is checked for its objective
 * value and its feasibility, and so are the valid bounds on the optimal
 * value. A sequence of random changes, with a fixed seed, is finally
 * checked against a copy of the data kept by the tester.
 *
 * The checks hold under NDEBUG as well, not being assert(). The exit code
 * is the number of failed checks.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFBlock.h"

#include "FakeSolver.h"

#include <algorithm>

#include <cmath>

#include <cstdio>

#include <fstream>

#include <functional>

#include <iostream>

#include <random>

#include <sstream>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using Range = Block::Range;
using Subset = Block::Subset;
using Vec_FNumber = MCFBlock::Vec_FNumber;
using Vec_CNumber = MCFBlock::Vec_CNumber;
using c_Vec_FNumber = MCFBlock::c_Vec_FNumber;
using c_Vec_CNumber = MCFBlock::c_Vec_CNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- GLOBALS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static int failures = 0;  // number of failed checks

static const double INF = Inf< double >();

/* The instance: 4 nodes, 4 units of flow from node 1 to node 4 (costs,
 * capacities):
 *
 *   0: 1 -> 2  (1, 3)    2: 2 -> 3  (1, 2)    4: 3 -> 4  (1, 5)
 *   1: 1 -> 3  (4, 5)    3: 2 -> 4  (5, 3)
 *
 * The optimum sends 2 units along 1-2-3-4 and 2 along 1-3-4, i.e., 16. */

static const Subset SN0 = { 1 , 1 , 2 , 2 , 3 };
static const Subset EN0 = { 2 , 3 , 3 , 4 , 4 };
static const Vec_CNumber C0 = { 1 , 4 , 1 , 5 , 1 };
static const Vec_FNumber U0 = { 3 , 5 , 2 , 3 , 5 };
static const Vec_FNumber B0 = { -4 , 0 , 0 , 4 };

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 }

/*--------------------------------------------------------------------------*/

static std::string str( double v )
{
 std::ostringstream s;
 s << v;
 return( s.str() );
 }

/*--------------------------------------------------------------------------*/
/// equality that takes two NaN as equal

static bool same( double a , double b )
{
 return( ( a == b ) || ( std::isnan( a ) && std::isnan( b ) ) );
 }

/*--------------------------------------------------------------------------*/
/// the instance above, dm of its arcs dynamic and room for mdm of them

static void load( MCFBlock & mcf , Index dm = 0 , Index mdm = 0 ,
		  c_Vec_FNumber & u = U0 , c_Vec_CNumber & c = C0 ,
		  c_Vec_FNumber & b = B0 , Index dn = 0 , Index mdn = 0 )
{
 mcf.load( 4 , 5 , EN0 , SN0 , u , c , b , dn , dm , mdn , mdm );
 }

/*--------------------------------------------------------------------------*/

static void generate( MCFBlock & mcf )
{
 mcf.generate_abstract_variables();
 mcf.generate_abstract_constraints();
 mcf.generate_objective();
 }

/*--------------------------------------------------------------------------*/

static LinearFunction * lf_of( FRowConstraint * cns )
{
 return( static_cast< LinearFunction * >( cns->get_function() ) );
 }

/*--------------------------------------------------------------------------*/
/// the abstract representation against the physical one
/** The flow Variable of every arc (fixed to 0 if the arc is closed or
 * deleted), the coefficients of the Objective (but that of a deleted arc,
 * whose flow is fixed to 0), the rhs of the bound Constraint and both sides
 * of the flow conservation ones, and the incidence of every arc, a deleted
 * one included, in the flow conservation Constraint of its two nodes. */

static void check_abstract( MCFBlock & mcf , const std::string & what ,
			    bool cons = true , bool bnds = true ,
			    bool obj = true )
{
 const Index m = mcf.get_NArcs();
 const Index n = mcf.get_NNodes();

 for( Index a = 0 ; a < m ; ++a ) {
  auto xa = mcf.i2p_x( a );
  check( mcf.p2i_x( xa ) == a , what + ": p2i_x( i2p_x( " +
	 std::to_string( a ) + " ) )" );
  if( mcf.is_deleted( a ) || mcf.is_closed( a ) )
   check( xa->is_fixed() && ( xa->get_value() == 0 ) , what + ": arc " +
	  std::to_string( a ) + " closed or deleted, flow not fixed to 0" );
  else
   check( ! xa->is_fixed() , what + ": arc " + std::to_string( a ) +
	  " open, flow fixed" );
  }

 if( obj ) {
  auto lf = static_cast< LinearFunction * >(
			  mcf.get_objective< FRealObjective >()->get_function() );
  check( lf->get_num_active_var() == m , what + ": " +
	 std::to_string( lf->get_num_active_var() ) +
	 " Variable in the Objective, " + std::to_string( m ) + " arcs" );
  for( Index a = 0 ; a < m ; ++a ) {
   const auto i = lf->is_active( mcf.i2p_x( a ) );
   check( i == a , what + ": arc " + std::to_string( a ) +
	  " in position " + std::to_string( i ) + " of the Objective" );
   // a deleted arc keeps whatever coefficient, its flow being fixed to 0
   if( ( i >= lf->get_num_active_var() ) || mcf.is_deleted( a ) )
    continue;
   const double ca = mcf.get_C( a );
   check( lf->get_coefficient( i ) == ca , what + ": arc " +
	  std::to_string( a ) + " has cost " + str( ca ) + " and coefficient "
	  + str( lf->get_coefficient( i ) ) );
   }
  }

 if( ! cons )
  return;

 if( bnds )
  for( Index a = 0 ; a < m ; ++a ) {
   if( mcf.is_deleted( a ) )
    continue;
   auto ub = mcf.i2p_ub( a );
   check( ub->get_active_var( 0 ) == mcf.i2p_x( a ) , what +
	  ": bound Constraint of arc " + std::to_string( a ) );
   check( ub->get_rhs() == mcf.get_U( a ) , what + ": arc " +
	  std::to_string( a ) + " has capacity " + str( mcf.get_U( a ) ) +
	  " and bound " + str( ub->get_rhs() ) );
   }

 Subset deg( n , 0 );
 for( Index a = 0 ; a < m ; ++a ) {
  ++deg[ mcf.get_SN( a ) - 1 ];
  ++deg[ mcf.get_EN( a ) - 1 ];
  }

 for( Index i = 0 ; i < n ; ++i ) {
  auto e = mcf.i2p_e( i );
  check( ( e->get_lhs() == mcf.get_B( i ) ) &&
	 ( e->get_rhs() == mcf.get_B( i ) ) , what + ": node " +
	 std::to_string( i ) + " has deficit " + str( mcf.get_B( i ) ) +
	 " and sides " + str( e->get_lhs() ) + ", " + str( e->get_rhs() ) );
  auto lf = lf_of( e );
  check( lf->get_num_active_var() == deg[ i ] , what + ": node " +
	 std::to_string( i ) + " has " + std::to_string( deg[ i ] ) +
	 " arcs and " + std::to_string( lf->get_num_active_var() ) +
	 " Variable in its flow conservation Constraint" );
  for( Index k = 0 ; k < lf->get_num_active_var() ; ++k ) {
   const auto a = mcf.p2i_x( lf->get_active_var( k ) );
   const auto cf = lf->get_coefficient( k );
   check( ( ( mcf.get_SN( a ) == i + 1 ) && ( cf == -1 ) ) ||
	  ( ( mcf.get_EN( a ) == i + 1 ) && ( cf == 1 ) ) , what +
	  ": arc " + std::to_string( a ) + " with coefficient " + str( cf ) +
	  " in the flow conservation Constraint of node " +
	  std::to_string( i ) );
   }
  }
 }  // end( check_abstract )

/*--------------------------------------------------------------------------*/
/*------------------------- THE MODIFICATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Modification received by a FakeSolver, GroupModification unpacked

struct Mods {
 std::vector< std::shared_ptr< MCFBlockMod > > phys;  ///< the physical ones
 Vec_sp_Mod abst;                                     ///< all the others
 };

static void flatten( const sp_Mod & mod , Mods & md )
{
 if( auto gm = std::dynamic_pointer_cast< GroupModification >( mod ) ) {
  for( auto & sm : gm->sub_Modifications() )
   flatten( sm , md );
  return;
  }

 if( auto pm = std::dynamic_pointer_cast< MCFBlockMod >( mod ) )
  md.phys.push_back( pm );
 else
  md.abst.push_back( mod );
 }

/*--------------------------------------------------------------------------*/
/// takes (and clears) the Modification the FakeSolver has received so far

static Mods take( FakeSolver * fs )
{
 Mods md;
 fs->lock_Modification_list();
 for( auto & mod : fs->get_Modification_list() )
  flatten( mod , md );
 fs->get_Modification_list().clear();
 fs->unlock_Modification_list();
 return( md );
 }

/*--------------------------------------------------------------------------*/
/// the number of abstract Modification of class T

template< class T >
static Index count( const Mods & md )
{
 Index k = 0;
 for( auto & mod : md.abst )
  if( std::dynamic_pointer_cast< T >( mod ) )
   ++k;
 return( k );
 }

/*--------------------------------------------------------------------------*/
/// no abstract Modification received still asks the Block to handle it

static void check_no_concern( const Mods & md , const std::string & what )
{
 for( auto & mod : md.abst )
  check( ! mod->concerns_Block() , what +
	 ": an abstract Modification received with concerns_Block()" );
 }

/*--------------------------------------------------------------------------*/
/// exactly one physical Modification, of type t, on the Range r

static void check_rngd( const Mods & md , int t , Range r ,
			const std::string & what )
{
 check( md.phys.size() == 1 , what + ": " + std::to_string( md.phys.size() )
	+ " physical Modification instead of 1" );
 if( md.phys.size() != 1 )
  return;
 auto rm = std::dynamic_pointer_cast< MCFBlockRngdMod >( md.phys.front() );
 check( rm && ( rm->type() == t ) && ( rm->rng() == r ) , what +
	": the physical Modification is not a MCFBlockRngdMod of type " +
	std::to_string( t ) + " on [ " + std::to_string( r.first ) + " , " +
	std::to_string( r.second ) + " )" );
 check_no_concern( md , what );
 }

/*--------------------------------------------------------------------------*/
/// exactly one physical Modification, of type t, on the Subset s

static void check_sbst( const Mods & md , int t , const Subset & s ,
			const std::string & what )
{
 check( md.phys.size() == 1 , what + ": " + std::to_string( md.phys.size() )
	+ " physical Modification instead of 1" );
 if( md.phys.size() != 1 )
  return;
 auto sm = std::dynamic_pointer_cast< MCFBlockSbstMod >( md.phys.front() );
 check( sm && ( sm->type() == t ) && ( sm->nms() == s ) , what +
	": the physical Modification is not a MCFBlockSbstMod of type " +
	std::to_string( t ) + " on the expected ordered Subset" );
 check_no_concern( md , what );
 }

/*--------------------------------------------------------------------------*/
/// nothing received at all

static void check_none( const Mods & md , const std::string & what )
{
 check( md.phys.empty() && md.abst.empty() , what + ": " +
	std::to_string( md.phys.size() ) + " physical and " +
	std::to_string( md.abst.size() ) +
	" abstract Modification where none was expected" );
 }

/*--------------------------------------------------------------------------*/
/// the MCFBlock with everything generated and a FakeSolver registered

static FakeSolver * attach( MCFBlock & mcf )
{
 generate( mcf );
 auto fs = new FakeSolver();
 mcf.register_Solver( fs );
 take( fs );
 return( fs );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- TESTS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/// load() with no capacities and no deficits: infinite and zero

static void test_load_no_data( void )
{
 MCFBlock mcf;
 // two dynamic nodes and two dynamic arcs yet to come
 mcf.load( 4 , 5 , EN0 , SN0 , {} , {} , {} , 0 , 0 , 2 , 2 );
 check( ( mcf.get_MaxNNodes() == 6 ) && ( mcf.get_MaxNArcs() == 7 ) ,
	"no data: maximum sizes" );
 check( mcf.get_U().empty() && mcf.get_B().empty() , "no data: U and B" );
 for( Index a = 0 ; a < 5 ; ++a )
  check( ( mcf.get_U( a ) == INF ) && ( mcf.get_C( a ) == 0 ) ,
	 "no data: capacity and cost of arc " + std::to_string( a ) );
 for( Index i = 0 ; i < 4 ; ++i )
  check( mcf.get_B( i ) == 0 , "no data: deficit of node " +
	 std::to_string( i ) );

 generate( mcf );
 check_abstract( mcf , "no data" );

 mcf.chg_ucap( INF , 0 );  // nothing to change
 check( mcf.get_U().empty() , "no data: infinite capacity makes U" );

 mcf.chg_ucap( 2 , 1 );
 check( mcf.get_U().size() == mcf.get_MaxNArcs() , "no data: chg_ucap() "
	"sizes U on " + std::to_string( mcf.get_U().size() ) );
 check( ( mcf.get_U( 1 ) == 2 ) && ( mcf.get_U( 0 ) == INF ) ,
	"no data: chg_ucap() values" );

 mcf.chg_dfct( 0 , 0 );  // nothing to change
 check( mcf.get_B().empty() , "no data: zero deficit makes B" );

 mcf.chg_dfct( -3 , 0 );
 mcf.chg_dfct( 3 , 3 );
 check( mcf.get_B().size() == mcf.get_MaxNNodes() , "no data: chg_dfct() "
	"sizes B on " + std::to_string( mcf.get_B().size() ) );
 check( ( mcf.get_B( 0 ) == -3 ) && ( mcf.get_B( 3 ) == 3 ) &&
	( mcf.get_B( 1 ) == 0 ) , "no data: chg_dfct() values" );
 check_abstract( mcf , "no data, changed" );

 // the vectors given only with infinite capacities and zero deficits
 MCFBlock mcf2;
 load( mcf2 , 0 , 0 , Vec_FNumber( 5 , INF ) , C0 , Vec_FNumber( 4 , 0 ) );
 check( mcf2.get_U().empty() && mcf2.get_B().empty() ,
	"all infinite, all zero: U and B" );
 }

/*--------------------------------------------------------------------------*/
/// the flow Variable of an arc of infinite cost is fixed, its cost zeroed

static void test_load_infinite_cost( void )
{
 MCFBlock mcf;
 load( mcf , 0 , 0 , U0 , { 1 , INF , 1 , 5 , 1 } );
 check( mcf.is_closed( 1 ) && ( mcf.get_C( 1 ) == 0 ) ,
	"infinite cost: arc closed, cost zeroed" );
 generate( mcf );
 check_abstract( mcf , "infinite cost" );
 }

/*--------------------------------------------------------------------------*/
/// load() refuses what it documents as wrong

static void test_load_errors( void )
{
 MCFBlock mcf;
 auto throws = [ & ]( const std::function< void( void ) > & f ) {
  try { f(); } catch( std::invalid_argument & ) { return( true ); }
  return( false );
  };

 check( throws( [ & ]() { mcf.load( 4 , 5 , EN0 , Subset( 4 , 1 ) ); } ) ,
	"load(): short pSn" );
 check( throws( [ & ]() { mcf.load( 4 , 5 , EN0 , SN0 , { 1 } ); } ) ,
	"load(): short pU" );
 check( throws( [ & ]() { mcf.load( 4 , 5 , EN0 , SN0 , {} , { 1 } ); } ) ,
	"load(): short pC" );
 check( throws( [ & ]() { mcf.load( 4 , 5 , EN0 , SN0 , {} , {} , { 1 } ); }
		) , "load(): short pB" );
 check( throws( [ & ]() { mcf.load( 3 , 5 , EN0 , SN0 ); } ) ,
	"load(): node out of range" );

 std::istringstream dmx( "p min 2 1\na 1 1 0 1 1\n" );
 check( throws( [ & ]() { mcf.load( dmx ); } ) , "load( DMX ): self-loop" );
 }

/*--------------------------------------------------------------------------*/
/// the round-trip instance: unbalanced deficits, an infinite capacity,
/// arcs 3 and 4 dynamic with room for one more, arc 0 closed, arc 3 deleted

static void load_round_trip( MCFBlock & mcf )
{
 load( mcf , 2 , 3 , { 3 , INF , 2 , 3 , 5 } , C0 , { -5 , 0 , 0 , 4 } );
 mcf.close_arc( 0 );
 mcf.remove_arc( 3 );
 }

/*--------------------------------------------------------------------------*/
/// print( 'C' ) and load( std::istream ) give the same instance back

static void test_dmx_round_trip( void )
{
 const std::string file = "MCFBlock_unit_test.dmx";

 MCFBlock mcf;
 load_round_trip( mcf );
 {
  std::ofstream out( file );
  mcf.print( out , 'C' );
  }

 MCFBlock mcf2;
 try {
  std::ifstream in( file );
  mcf2.load( in );
  }
 catch( std::exception & e ) {
  check( false , std::string( "DMX: load() of what print() wrote: " ) +
	 e.what() );
  std::remove( file.c_str() );
  return;
  }
 std::remove( file.c_str() );

 // a DMX file is fully static, and has no closed or deleted arcs: those are
 // written as arcs of zero capacity
 check( ( mcf2.get_NNodes() == 4 ) && ( mcf2.get_NArcs() == 5 ) &&
	( mcf2.get_NStaticArcs() == 5 ) , "DMX: sizes" );
 const Vec_CNumber C = { 1 , 4 , 1 , 1e6 , 1 };
 const Vec_FNumber U = { 0 , INF , 2 , 0 , 5 };
 for( Index a = 0 ; a < mcf2.get_NArcs() ; ++a ) {
  check( ( mcf2.get_SN( a ) == SN0[ a ] ) && ( mcf2.get_EN( a ) == EN0[ a ] )
	 , "DMX: nodes of arc " + std::to_string( a ) );
  check( ! mcf2.is_closed( a ) , "DMX: arc " + std::to_string( a ) +
	 " closed after the round trip" );
  check( mcf2.get_C( a ) == C[ a ] , "DMX: cost of arc " +
	 std::to_string( a ) + " is " + str( mcf2.get_C( a ) ) );
  check( mcf2.get_U( a ) == U[ a ] , "DMX: capacity of arc " +
	 std::to_string( a ) + " is " + str( mcf2.get_U( a ) ) );
  }
 for( Index i = 0 ; i < 4 ; ++i )
  check( mcf2.get_B( i ) == mcf.get_B( i ) , "DMX: deficit of node " +
	 std::to_string( i ) );

 // the capacities given as INF, Inf, inf and +Inf
 std::istringstream dmx( "c a comment\np min 3 4\nn 1 2\nn 3 -2\n"
			 "a 1 2 0 INF 1\na 2 3 0 Inf 1\na 1 3 0 inf 3\n"
			 "a 1 3 1 +Inf 3\n" );
 MCFBlock mcf3;
 try {
  mcf3.load( dmx );
  check( mcf3.get_U().empty() , "DMX: infinite capacities" );
  // the lower bound of 1 moves one unit from node 1 to node 3
  check( ( mcf3.get_B( 0 ) == -1 ) && ( mcf3.get_B( 2 ) == 1 ) ,
	 "DMX: deficits and lower bound" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "DMX: infinite capacities: " ) + e.what() );
  }
 }

/*--------------------------------------------------------------------------*/
/// serialize() and deserialize() give the same instance back

static void test_netcdf_round_trip( void )
{
 const std::string file = "MCFBlock_unit_test.nc4";

 MCFBlock mcf;
 load_round_trip( mcf );
 {
  netCDF::NcFile f( file , netCDF::NcFile::replace );
  mcf.serialize( f );
  }

 MCFBlock mcf2;
 {
  netCDF::NcFile f( file , netCDF::NcFile::read );
  mcf2.deserialize( f );
  }
 std::remove( file.c_str() );

 check( ( mcf2.get_NNodes() == 4 ) && ( mcf2.get_NArcs() == 5 ) &&
	( mcf2.get_NStaticArcs() == 3 ) && ( mcf2.get_MaxNArcs() == 6 ) &&
	( mcf2.get_MaxNNodes() == 4 ) , "netCDF: sizes" );
 for( Index a = 0 ; a < mcf2.get_NArcs() ; ++a ) {
  check( ( mcf2.get_SN( a ) == SN0[ a ] ) && ( mcf2.get_EN( a ) == EN0[ a ] )
	 , "netCDF: nodes of arc " + std::to_string( a ) );
  check( mcf2.is_deleted( a ) == mcf.is_deleted( a ) , "netCDF: deletion "
	 "of arc " + std::to_string( a ) );
  // a closed arc comes back open with capacity 0 and the same cost, as in
  // the DMX file [see test_dmx_round_trip()]
  check( ! mcf2.is_closed( a ) , "netCDF: arc " + std::to_string( a ) +
	 " closed after the round trip" );
  if( mcf.is_deleted( a ) )
   continue;
  check( mcf2.get_C( a ) == mcf.get_C( a ) , "netCDF: cost of arc " +
	 std::to_string( a ) );
  const double ua = mcf.is_closed( a ) ? 0 : mcf.get_U( a );
  check( mcf2.get_U( a ) == ua , "netCDF: capacity of arc " +
	 std::to_string( a ) + " is " + str( mcf2.get_U( a ) ) +
	 " instead of " + str( ua ) );
  }
 for( Index i = 0 ; i < 4 ; ++i )
  check( mcf2.get_B( i ) == mcf.get_B( i ) , "netCDF: deficit of node " +
	 std::to_string( i ) );

 generate( mcf2 );
 check_abstract( mcf2 , "netCDF" );

 // the deleted arc takes the next new arc, as it did before
 check( mcf2.add_arc( 1 , 4 , 2 , 7 ) == 3 , "netCDF: name of the arc added "
	"in the deleted slot" );
 check_abstract( mcf2 , "netCDF, arc added" );

 // no capacities and a closed arc: U is written for the zero, in both
 // formats, which agree
 MCFBlock mcf5;
 load( mcf5 , 0 , 0 , {} );
 mcf5.close_arc( 2 );
 {
  netCDF::NcFile f( file , netCDF::NcFile::replace );
  mcf5.serialize( f );
  }
 MCFBlock mcf6;
 {
  netCDF::NcFile f( file , netCDF::NcFile::read );
  mcf6.deserialize( f );
  }
 std::remove( file.c_str() );
 MCFBlock mcf7;
 try {
  std::ostringstream out;
  mcf5.print( out , 'C' );
  std::istringstream in( out.str() );
  mcf7.load( in );
  }
 catch( std::exception & e ) {
  check( false , std::string( "no capacities, closed arc: DMX: " ) +
	 e.what() );
  }
 for( Index a = 0 ; a < 5 ; ++a ) {
  const double ua = a == 2 ? 0 : INF;
  check( ( mcf6.get_U( a ) == ua ) && ( ! mcf6.is_closed( a ) ) &&
	 ( mcf6.get_C( a ) == C0[ a ] ) , "no capacities, closed arc: netCDF: "
	 "arc " + std::to_string( a ) + " has capacity " +
	 str( mcf6.get_U( a ) ) );
  check( ( mcf7.get_U( a ) == mcf6.get_U( a ) ) &&
	 ( mcf7.get_C( a ) == mcf6.get_C( a ) ) &&
	 ( mcf7.get_SN( a ) == mcf6.get_SN( a ) ) &&
	 ( mcf7.get_EN( a ) == mcf6.get_EN( a ) ) , "no capacities, closed "
	 "arc: DMX and netCDF disagree on arc " + std::to_string( a ) );
  }
 for( Index i = 0 ; i < 4 ; ++i )
  check( ( mcf6.get_B( i ) == B0[ i ] ) && ( mcf7.get_B( i ) == B0[ i ] ) ,
	 "no capacities, closed arc: deficit of node " + std::to_string( i ) );

 // no capacities and no deficits stay so
 MCFBlock mcf3;
 load( mcf3 , 0 , 0 , {} , C0 , {} );
 {
  netCDF::NcFile f( file , netCDF::NcFile::replace );
  mcf3.serialize( f );
  }
 MCFBlock mcf4;
 {
  netCDF::NcFile f( file , netCDF::NcFile::read );
  mcf4.deserialize( f );
  }
 std::remove( file.c_str() );
 check( mcf4.get_U().empty() && mcf4.get_B().empty() &&
	( mcf4.get_NArcs() == 5 ) , "netCDF: no capacities, no deficits" );
 }

/*--------------------------------------------------------------------------*/
/// an instance with no arcs through both formats

static void test_zero_arcs( void )
{
 MCFBlock mcf;
 mcf.load( 3 , 0 , {} , {} , {} , {} , { -1 , 0 , 1 } );
 generate( mcf );
 check_abstract( mcf , "zero arcs" );
 check( ! mcf.is_feasible() , "zero arcs: feasible with deficits" );

 std::ostringstream out;
 mcf.print( out , 'C' );
 MCFBlock mcf2;
 try {
  std::istringstream in( out.str() );
  mcf2.load( in );
  check( ( mcf2.get_NNodes() == 3 ) && ( mcf2.get_NArcs() == 0 ) &&
	 ( mcf2.get_B( 0 ) == -1 ) && ( mcf2.get_B( 2 ) == 1 ) ,
	 "zero arcs: DMX" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "zero arcs: DMX: " ) + e.what() );
  }

 const std::string file = "MCFBlock_test_0.nc4";
 MCFBlock mcf3;
 try {
  {
   netCDF::NcFile f( file , netCDF::NcFile::replace );
   mcf.serialize( f );
   }
  netCDF::NcFile f( file , netCDF::NcFile::read );
  mcf3.deserialize( f );
  check( ( mcf3.get_NNodes() == 3 ) && ( mcf3.get_NArcs() == 0 ) &&
	 ( mcf3.get_B( 0 ) == -1 ) && ( mcf3.get_B( 2 ) == 1 ) ,
	 "zero arcs: netCDF" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "zero arcs: netCDF: " ) + e.what() );
  }
 std::remove( file.c_str() );
 }

/*--------------------------------------------------------------------------*/
/// the costs, in every form

static void test_costs( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );
 auto fs = attach( mcf );

 mcf.chg_cost( 10 , 2 );
 auto md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCost , Range( 2 , 3 ) , "chg_cost()" );
 check( count< C05FunctionModLin >( md ) == 1 , "chg_cost(): abstract" );
 check( mcf.get_C( 2 ) == 10 , "chg_cost(): value" );
 check_abstract( mcf , "chg_cost()" );

 mcf.chg_cost( 10 , 2 );  // no change
 check_none( take( fs ) , "chg_cost() to the same value" );

 const Vec_CNumber nc = { 7 , 8 };
 mcf.chg_costs( nc.begin() , Range( 0 , 2 ) );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCost , Range( 0 , 2 ) , "chg_costs( Range )"
	     );
 check( ( mcf.get_C( 0 ) == 7 ) && ( mcf.get_C( 1 ) == 8 ) ,
	"chg_costs( Range ): values" );
 check_abstract( mcf , "chg_costs( Range )" );

 // a Range beyond the arcs is cut, the costs of the rest are not read
 const Vec_CNumber nc2 = { 6 , 6 };
 mcf.chg_costs( std::span< const double >( nc2 ) , Range( 3 , 10 ) );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCost , Range( 3 , 5 ) ,
	     "chg_costs( Range beyond the arcs )" );
 check_abstract( mcf , "chg_costs( Range beyond the arcs )" );

 mcf.chg_costs( nc.begin() , Range( 3 , 3 ) );
 check_none( take( fs ) , "chg_costs( empty Range )" );
 mcf.chg_costs( nc.begin() , Subset() );
 check_none( take( fs ) , "chg_costs( empty Subset )" );

 // a short span is refused
 bool thrown = false;
 try {
  mcf.chg_costs( std::span< const double >( nc ) , Range( 0 , 3 ) );
  }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "chg_costs(): span shorter than the Range" );

 const Vec_CNumber nc3 = { 9 , 3 };
 mcf.chg_costs( nc3.begin() , Subset( { 4 , 1 } ) );  // unordered
 md = take( fs );
 check_sbst( md , MCFBlockMod::eChgCost , Subset( { 1 , 4 } ) ,
	     "chg_costs( Subset )" );
 check( ( mcf.get_C( 4 ) == 9 ) && ( mcf.get_C( 1 ) == 3 ) ,
	"chg_costs( Subset ): values" );
 check_abstract( mcf , "chg_costs( Subset )" );

 // the cost of a closed arc changes, and holds when it is opened
 mcf.close_arc( 0 );
 take( fs );
 mcf.chg_cost( 2 , 0 );
 mcf.open_arc( 0 );
 take( fs );
 check( mcf.get_C( 0 ) == 2 , "chg_cost() of a closed arc" );
 check_abstract( mcf , "chg_cost() of a closed arc" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// the capacities, in every form

static void test_capacities( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );
 auto fs = attach( mcf );

 mcf.chg_ucap( 4 , 3 );
 auto md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCaps , Range( 3 , 4 ) , "chg_ucap()" );
 check( count< RowConstraintMod >( md ) == 1 , "chg_ucap(): abstract" );
 check_abstract( mcf , "chg_ucap()" );

 const Vec_FNumber nu = { 6 , INF , 1 };
 mcf.chg_ucaps( nu.begin() , Range( 1 , 4 ) );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCaps , Range( 1 , 4 ) ,
	     "chg_ucaps( Range )" );
 check( count< RowConstraintMod >( md ) == 3 , "chg_ucaps( Range ): " +
	std::to_string( count< RowConstraintMod >( md ) ) + " abstract" );
 check( ( mcf.get_U( 1 ) == 6 ) && ( mcf.get_U( 2 ) == INF ) &&
	( mcf.get_U( 3 ) == 1 ) , "chg_ucaps( Range ): values" );
 check_abstract( mcf , "chg_ucaps( Range )" );

 mcf.chg_ucaps( nu.begin() , Range( 2 , 2 ) );
 check_none( take( fs ) , "chg_ucaps( empty Range )" );
 mcf.chg_ucaps( nu.begin() , Subset() );
 check_none( take( fs ) , "chg_ucaps( empty Subset )" );

 const Vec_FNumber nu2 = { 0 , 8 };
 mcf.chg_ucaps( nu2.begin() , Subset( { 4 , 0 } ) );  // unordered
 md = take( fs );
 check_sbst( md , MCFBlockMod::eChgCaps , Subset( { 0 , 4 } ) ,
	     "chg_ucaps( Subset )" );
 check( ( mcf.get_U( 4 ) == 0 ) && ( mcf.get_U( 0 ) == 8 ) ,
	"chg_ucaps( Subset ): values" );
 check_abstract( mcf , "chg_ucaps( Subset )" );

 mcf.unregister_Solvers( true );

 // without the bound Constraint a finite capacity cannot be set
 MCFBlock mcf2;
 load( mcf2 , 0 , 0 , {} );
 auto nobnd = new SimpleConfiguration< int >( 1 );
 mcf2.generate_abstract_constraints( nobnd );
 delete nobnd;
 bool thrown = false;
 try { mcf2.chg_ucap( 1 , 0 ); } catch( std::logic_error & ) { thrown = true; }
 check( thrown , "chg_ucap() without bound Constraint" );
 }

/*--------------------------------------------------------------------------*/
/// the deficits, in every form, dynamic nodes included

static void test_deficits( void )
{
 MCFBlock mcf;
 load( mcf , 0 , 0 , U0 , C0 , B0 , 1 );  // node 3 dynamic
 check( mcf.get_NStaticNodes() == 3 , "dynamic node" );
 auto fs = attach( mcf );

 mcf.chg_dfct( 1 , 3 );
 auto md = take( fs );
 check_rngd( md , MCFBlockMod::eChgDfct , Range( 3 , 4 ) ,
	     "chg_dfct() of the dynamic node" );
 check( count< RowConstraintMod >( md ) == 1 , "chg_dfct(): abstract" );
 check_abstract( mcf , "chg_dfct()" );

 const Vec_FNumber nb = { -2 , 1 , 1 };
 mcf.chg_dfcts( nb.begin() , Range( 1 , 4 ) );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eChgDfct , Range( 1 , 4 ) ,
	     "chg_dfcts( Range )" );
 check_abstract( mcf , "chg_dfcts( Range )" );

 mcf.chg_dfcts( nb.begin() , Range( 4 , 9 ) );
 check_none( take( fs ) , "chg_dfcts( Range beyond the nodes )" );
 mcf.chg_dfcts( nb.begin() , Subset() );
 check_none( take( fs ) , "chg_dfcts( empty Subset )" );

 // unordered, the dynamic node first: its name is beyond the static nodes
 // but not beyond the static arcs
 const Vec_FNumber nb2 = { 5 , -5 };
 mcf.chg_dfcts( nb2.begin() , Subset( { 3 , 0 } ) );
 md = take( fs );
 check_sbst( md , MCFBlockMod::eChgDfct , Subset( { 0 , 3 } ) ,
	     "chg_dfcts( unordered Subset )" );
 check( ( mcf.get_B( 0 ) == -5 ) && ( mcf.get_B( 3 ) == 5 ) ,
	"chg_dfcts( unordered Subset ): values" );
 check_abstract( mcf , "chg_dfcts( unordered Subset )" );

 const Vec_FNumber nb3 = { -4 , 4 };
 mcf.chg_dfcts( nb3.begin() , Subset( { 0 , 3 } ) , true );
 md = take( fs );
 check_sbst( md , MCFBlockMod::eChgDfct , Subset( { 0 , 3 } ) ,
	     "chg_dfcts( ordered Subset )" );
 check_abstract( mcf , "chg_dfcts( ordered Subset )" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// closing and opening arcs, in every form

static void test_close_open( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );
 auto fs = attach( mcf );

 mcf.close_arcs( Range( 1 , 4 ) );  // static and dynamic
 auto md = take( fs );
 check_rngd( md , MCFBlockMod::eCloseArc , Range( 1 , 4 ) ,
	     "close_arcs( Range )" );
 check( count< VariableMod >( md ) == 3 , "close_arcs( Range ): abstract" );
 for( Index a = 1 ; a < 4 ; ++a )
  check( mcf.is_closed( a ) , "close_arcs( Range ): arc " +
	 std::to_string( a ) );
 check_abstract( mcf , "close_arcs( Range )" );

 mcf.close_arcs( Range( 1 , 4 ) );  // closed already
 check_none( take( fs ) , "close_arcs() of closed arcs" );
 mcf.close_arcs( Range( 2 , 2 ) );
 check_none( take( fs ) , "close_arcs( empty Range )" );
 mcf.close_arcs( Subset() );
 check_none( take( fs ) , "close_arcs( empty Subset )" );

 mcf.open_arcs( Subset( { 3 , 1 } ) );  // unordered
 md = take( fs );
 check_sbst( md , MCFBlockMod::eOpenArc , Subset( { 1 , 3 } ) ,
	     "open_arcs( Subset )" );
 check( ( ! mcf.is_closed( 1 ) ) && mcf.is_closed( 2 ) &&
	( ! mcf.is_closed( 3 ) ) , "open_arcs( Subset ): arcs" );
 check_abstract( mcf , "open_arcs( Subset )" );

 mcf.open_arcs( Range( 0 , 2 ) );  // open already
 check_none( take( fs ) , "open_arcs() of open arcs" );
 mcf.open_arcs( Subset() );
 check_none( take( fs ) , "open_arcs( empty Subset )" );

 mcf.close_arcs( Subset( { 4 , 0 } ) );
 md = take( fs );
 check_sbst( md , MCFBlockMod::eCloseArc , Subset( { 0 , 4 } ) ,
	     "close_arcs( Subset )" );
 mcf.open_arcs( Range( 0 , 10 ) );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eOpenArc , Range( 0 , 5 ) ,
	     "open_arcs( Range beyond the arcs )" );
 check_abstract( mcf , "open_arcs( Range )" );

 mcf.close_arc( 4 );
 check_rngd( take( fs ) , MCFBlockMod::eCloseArc , Range( 4 , 5 ) ,
	     "close_arc()" );
 mcf.open_arc( 4 );
 check_rngd( take( fs ) , MCFBlockMod::eOpenArc , Range( 4 , 5 ) ,
	     "open_arc()" );
 check_abstract( mcf , "close_arc(), open_arc()" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// eNoMod does the change and issues nothing, eModBlck issues each
/// Modification once and none of them asks the MCFBlock to handle it

static void test_mod_params( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );
 auto fs = attach( mcf );

 // the changes are made, nothing is issued
 mcf.chg_cost( 3 , 0 , eNoMod , eNoMod );
 mcf.chg_ucap( 1 , 4 , eNoMod , eNoMod );
 mcf.chg_dfct( 1 , 1 , eNoMod , eNoMod );
 mcf.close_arc( 1 , eNoMod , eNoMod );
 check_none( take( fs ) , "eNoMod [Observer::issue_pmod() takes eNoMod for "
	     "a request to issue the physical Modification]" );
 check( ( mcf.get_C( 0 ) == 3 ) && ( mcf.get_U( 4 ) == 1 ) &&
	( mcf.get_B( 1 ) == 1 ) && mcf.is_closed( 1 ) , "eNoMod: values" );
 check_abstract( mcf , "eNoMod" );
 mcf.chg_dfct( 0 , 1 , eNoMod , eNoMod );
 mcf.open_arc( 1 , eNoMod , eNoMod );
 take( fs );

 mcf.chg_ucap( 4 , 0 , eModBlck , eModBlck );  // static arc
 auto md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCaps , Range( 0 , 1 ) ,
	     "chg_ucap( eModBlck ) of a static arc" );
 check( count< RowConstraintMod >( md ) == 1 , "chg_ucap( eModBlck ): " +
	std::to_string( count< RowConstraintMod >( md ) ) + " abstract" );

 mcf.chg_ucap( 4 , 3 , eModBlck , eModBlck );  // dynamic arc
 check_rngd( take( fs ) , MCFBlockMod::eChgCaps , Range( 3 , 4 ) ,
	     "chg_ucap( eModBlck ) of a dynamic arc" );

 mcf.chg_dfct( 2 , 1 , eModBlck , eModBlck );
 check_rngd( take( fs ) , MCFBlockMod::eChgDfct , Range( 1 , 2 ) ,
	     "chg_dfct( eModBlck )" );
 mcf.chg_cost( 2 , 1 , eModBlck , eModBlck );
 check_rngd( take( fs ) , MCFBlockMod::eChgCost , Range( 1 , 2 ) ,
	     "chg_cost( eModBlck )" );
 mcf.close_arc( 2 , eModBlck , eModBlck );
 check_rngd( take( fs ) , MCFBlockMod::eCloseArc , Range( 2 , 3 ) ,
	     "close_arc( eModBlck )" );
 check_abstract( mcf , "eModBlck" );

 // eDryRun for the abstract representation leaves it alone, the physical
 // Modification being issued all the same
 mcf.chg_cost( 5 , 4 , eNoBlck , eDryRun );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eChgCost , Range( 4 , 5 ) ,
	     "chg_cost( eDryRun )" );
 check( md.abst.empty() , "chg_cost( eDryRun ): abstract" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// adding and removing arcs

static void test_add_remove( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );  // arcs 3 and 4 dynamic, room for one more
 auto fs = attach( mcf );

 check( mcf.add_arc( 1 , 4 , 7 , 10 ) == 5 , "add_arc(): name" );
 auto md = take( fs );
 check_rngd( md , MCFBlockMod::eAddArc , Range( 5 , 6 ) , "add_arc()" );
 check( count< BlockModAD >( md ) >= 1 , "add_arc(): no BlockModAdd" );
 check( ( mcf.get_NArcs() == 6 ) && ( mcf.get_C( 5 ) == 7 ) &&
	( mcf.get_U( 5 ) == 10 ) , "add_arc(): data" );
 check_abstract( mcf , "add_arc()" );

 check( mcf.add_arc( 1 , 2 ) == Inf< Index >() , "add_arc() when full" );
 check_none( take( fs ) , "add_arc() when full" );

 bool thrown = false;
 try { mcf.remove_arc( 0 ); }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "remove_arc() of a static arc" );
 thrown = false;
 try { mcf.add_arc( 0 , 2 ); }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "add_arc() from node 0" );

 mcf.remove_arc( 3 );  // in the middle
 md = take( fs );
 check_rngd( md , MCFBlockMod::eRmvArc , Range( 3 , 4 ) ,
	     "remove_arc() in the middle" );
 check( mcf.is_deleted( 3 ) && ( mcf.get_NArcs() == 6 ) &&
	( ! mcf.is_closed( 3 ) ) , "remove_arc() in the middle: data" );
 check_abstract( mcf , "remove_arc() in the middle" );

 mcf.remove_arc( 3 );  // deleted already
 check_none( take( fs ) , "remove_arc() of a deleted arc" );
 mcf.chg_cost( 1 , 3 );
 check( mcf.is_deleted( 3 ) , "chg_cost() of a deleted arc" );
 mcf.open_arc( 3 );
 check( mcf.is_deleted( 3 ) , "open_arc() of a deleted arc" );
 take( fs );

 // the hole is filled, from other nodes
 check( mcf.add_arc( 1 , 4 , 2 , 1 ) == 3 , "add_arc() in the hole" );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eAddArc , Range( 3 , 4 ) ,
	     "add_arc() in the hole" );
 check( ( mcf.get_SN( 3 ) == 1 ) && ( mcf.get_EN( 3 ) == 4 ) &&
	( mcf.get_C( 3 ) == 2 ) && ( mcf.get_U( 3 ) == 1 ) ,
	"add_arc() in the hole: data" );
 check_abstract( mcf , "add_arc() in the hole" );

 // the last arc goes, and the deleted ones before it with it
 mcf.remove_arc( 4 );
 take( fs );
 mcf.remove_arc( 5 );
 md = take( fs );
 check_rngd( md , MCFBlockMod::eRmvArc , Range( 4 , 6 ) ,
	     "remove_arc() of the last arc" );
 check( mcf.get_NArcs() == 4 , "remove_arc() of the last arc: " +
	std::to_string( mcf.get_NArcs() ) + " arcs" );
 check_abstract( mcf , "remove_arc() of the last arc" );

 mcf.remove_arc( 3 );
 take( fs );
 check( mcf.get_NArcs() == 3 , "all dynamic arcs removed" );
 check_abstract( mcf , "all dynamic arcs removed" );

 check( mcf.add_arc( 3 , 1 , 1 , INF ) == 3 , "add_arc() after removing" );
 take( fs );
 check_abstract( mcf , "add_arc() after removing" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// arcs deleted before the flow conservation Constraint are generated

static void test_remove_before_generate( void )
{
 MCFBlock mcf;
 load( mcf , 3 , 3 );  // arcs 2, 3 and 4 dynamic
 mcf.remove_arc( 2 );
 mcf.remove_arc( 3 );
 generate( mcf );
 check_abstract( mcf , "deleted before generate" );

 try {
  check( mcf.add_arc( 1 , 4 , 3 , 2 ) == 2 , "deleted before generate: "
	 "add_arc() name" );
  check_abstract( mcf , "deleted before generate, add_arc() in the hole" );
  mcf.remove_arc( 4 );
  check( mcf.get_NArcs() == 3 , "deleted before generate: remove_arc() of "
	 "the last arc" );
  check_abstract( mcf , "deleted before generate, remove_arc()" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "deleted before generate: " ) + e.what() );
  }
 }

/*--------------------------------------------------------------------------*/
/// every arc dynamic, all of them removed

static void test_remove_all( void )
{
 MCFBlock mcf;
 load( mcf , 5 , 5 );
 check( mcf.get_NStaticArcs() == 0 , "all dynamic: static arcs" );
 generate( mcf );
 for( Index a = 0 ; a < 4 ; ++a )  // in the middle
  mcf.remove_arc( a );
 mcf.remove_arc( 4 );             // the last, and all the others with it
 check( mcf.get_NArcs() == 0 , "all dynamic, all removed: " +
	std::to_string( mcf.get_NArcs() ) + " arcs" );
 check_abstract( mcf , "all dynamic, all removed" );

 check( mcf.add_arc( 1 , 4 , 1 , 4 ) == 0 , "all removed: add_arc()" );
 check_abstract( mcf , "all removed, add_arc()" );
 mcf.remove_arc( 0 );
 check( mcf.get_NArcs() == 0 , "the only arc removed" );
 }

/*--------------------------------------------------------------------------*/
/// a Range of costs and capacities ending on a deleted arc

static void test_range_ending_deleted( void )
{
 MCFBlock mcf;
 load( mcf , 3 , 3 );  // arcs 2, 3 and 4 dynamic
 mcf.remove_arc( 3 );
 auto fs = attach( mcf );

 // arc 2 changes to the value arc 3 would have had, arc 3 is deleted
 const Vec_CNumber nc = { 9 , 2 , 1 };
 mcf.chg_costs( nc.begin() , Range( 1 , 4 ) );
 take( fs );
 check( ( mcf.get_C( 1 ) == 9 ) && ( mcf.get_C( 2 ) == 2 ) , "chg_costs() "
	"of a Range ending on a deleted arc: costs " + str( mcf.get_C( 1 ) ) +
	", " + str( mcf.get_C( 2 ) ) );
 check_abstract( mcf , "chg_costs() of a Range ending on a deleted arc" );

 const Vec_FNumber nu = { 9 , 5 , 2 };
 mcf.chg_ucaps( nu.begin() , Range( 1 , 4 ) );
 take( fs );
 check( ( mcf.get_U( 1 ) == 9 ) && ( mcf.get_U( 2 ) == 5 ) , "chg_ucaps() "
	"of a Range ending on a deleted arc: capacities " +
	str( mcf.get_U( 1 ) ) + ", " + str( mcf.get_U( 2 ) ) );
 check_abstract( mcf , "chg_ucaps() of a Range ending on a deleted arc" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// changing the abstract representation changes the data

static void test_abstract_to_physical( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );
 auto fs = attach( mcf );

 auto lfo = static_cast< LinearFunction * >(
			  mcf.get_objective< FRealObjective >()->get_function() );
 lfo->modify_coefficient( 2 , 11 );
 auto md = take( fs );
 check( mcf.get_C( 2 ) == 11 , "objective coefficient: cost" );
 check_rngd( md , MCFBlockMod::eChgCost , Range( 2 , 3 ) ,
	     "objective coefficient" );

 lfo->modify_coefficients( { 6 , 6 } , Range( 3 , 5 ) );
 take( fs );
 check( ( mcf.get_C( 3 ) == 6 ) && ( mcf.get_C( 4 ) == 6 ) ,
	"objective coefficients: costs" );

 mcf.i2p_ub( 4 )->set_rhs( 1 );
 md = take( fs );
 check( mcf.get_U( 4 ) == 1 , "bound rhs: capacity" );
 check_rngd( md , MCFBlockMod::eChgCaps , Range( 4 , 5 ) , "bound rhs" );

 mcf.i2p_e( 1 )->set_both( 2 );
 md = take( fs );
 check( mcf.get_B( 1 ) == 2 , "flow conservation sides: deficit" );
 check_rngd( md , MCFBlockMod::eChgDfct , Range( 1 , 2 ) ,
	     "flow conservation sides" );

 mcf.i2p_x( 0 )->set_value( 0 );
 mcf.i2p_x( 0 )->is_fixed( true );
 md = take( fs );
 check( mcf.is_closed( 0 ) , "fixing the flow: closed" );
 check_rngd( md , MCFBlockMod::eCloseArc , Range( 0 , 1 ) , "fixing the flow"
	     );
 mcf.i2p_x( 0 )->is_fixed( false );
 md = take( fs );
 check( ! mcf.is_closed( 0 ) , "unfixing the flow: open" );
 check_rngd( md , MCFBlockMod::eOpenArc , Range( 0 , 1 ) ,
	     "unfixing the flow" );
 check_abstract( mcf , "abstract to physical" );

 // what cannot be done through the abstract representation
 bool thrown = false;
 try { mcf.i2p_e( 1 )->set_rhs( 3 ); }
 catch( std::exception & ) { thrown = true; }
 check( thrown , "one side of a flow conservation Constraint" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// a flow set by hand: objective value and feasibility

static void test_flow_by_hand( void )
{
 MCFBlock mcf;
 load( mcf , 2 , 3 );
 generate( mcf );

 const Vec_FNumber opt = { 2 , 2 , 2 , 0 , 4 };
 mcf.set_x( opt.begin() );
 auto obj = mcf.get_objective< FRealObjective >();
 obj->compute();
 check( mcf.get_objective_value() == 16 , "objective value " +
	str( mcf.get_objective_value() ) );
 check( mcf.is_feasible() && mcf.is_feasible( true ) , "optimal flow" );

 Vec_FNumber got( 5 );
 mcf.get_x( got.begin() );
 check( got == opt , "get_x()" );
 mcf.get_x( got.begin() , Range( 4 , 5 ) );  // the last dynamic arc
 check( got[ 0 ] == 4 , "get_x( Range ) from the second dynamic arc gives "
	+ str( got[ 0 ] ) );
 mcf.get_x( got.begin() , Subset( { 1 , 4 } ) );
 check( ( got[ 0 ] == 2 ) && ( got[ 1 ] == 4 ) , "get_x( Subset )" );

 auto sol = dynamic_cast< MCFSolution * >( mcf.get_Solution( nullptr ,
							      false ) );
 check( sol && mcf.is_sol_feasible( sol ) , "is_sol_feasible()" );

 // conservation broken at nodes 3 and 4
 mcf.set_x( 4 , 3 );
 check( ( ! mcf.flow_feasible( 0 ) ) && mcf.bound_feasible( 0 ) &&
	( ! mcf.is_feasible() ) , "conservation broken" );
 check( ! mcf.flow_feasible( 0 , true ) , "conservation broken, abstract" );
 // off by 1 at node 3, of deficit 0, and by 1 / 4 at node 4, of deficit 4
 auto eps = new SimpleConfiguration< double >( 0.5 );
 check( ! mcf.is_feasible( false , eps ) , "conservation broken beyond the "
	"tolerance" );
 eps->f_value = 1.5;
 check( mcf.is_feasible( false , eps ) , "conservation broken within the "
	"tolerance" );
 delete eps;

 // 3 units along 1-2-3-4 and 1 along 1-3-4: arc 2 exceeds its capacity
 const Vec_FNumber over = { 3 , 1 , 3 , 0 , 4 };
 mcf.set_x( over.begin() );
 check( mcf.flow_feasible( 0 ) && ( ! mcf.bound_feasible( 0 ) ) ,
	"capacity exceeded" );
 check( ! mcf.bound_feasible( 0 , true ) , "capacity exceeded, abstract" );
 check( mcf.is_sol_feasible( sol ) , "the MCFSolution keeps its own flow" );

 // flow on a closed arc, which the Variable cannot tell
 check( mcf.bound_feasible( 0 , Vec_FNumber( { 2 , 2 , 2 , 0 , 4 } ) ) ,
	"open arcs, explicit flow" );
 mcf.close_arc( 3 );
 check( ! mcf.bound_feasible( 0 , Vec_FNumber( { 1 , 3 , 1 , 1 , 3 } ) ) ,
	"flow on a closed arc" );
 check( mcf.bound_feasible( 0 , Vec_FNumber( { 2 , 2 , 2 , 0 , 4 } ) ) ,
	"no flow on a closed arc" );

 // unbalanced deficits: no flow conserves
 mcf.chg_dfct( -5 , 0 );
 mcf.set_x( opt.begin() );
 check( ! mcf.flow_feasible( 0 ) , "unbalanced deficits" );

 // a deleted arc does not count
 mcf.chg_dfct( -4 , 0 );
 mcf.open_arc( 3 );
 mcf.remove_arc( 3 );
 check( mcf.flow_feasible( 0 , Vec_FNumber( { 2 , 2 , 2 , 9 , 4 } ) ) ,
	"flow on a deleted arc" );

 delete sol;
 }

/*--------------------------------------------------------------------------*/
/// the valid bounds on the optimal value

static void test_bounds( void )
{
 MCFBlock mcf;
 load( mcf );

 // sum of c * u over the positive costs, over the negative ones
 check( mcf.get_valid_upper_bound() == INF , "global upper bound" );
 check( mcf.get_valid_upper_bound( true ) == 45 , "conditional upper bound "
	"right after load(): " + str( mcf.get_valid_upper_bound( true ) ) );
 check( mcf.get_valid_lower_bound() == 0 , "lower bound after load(): " +
	str( mcf.get_valid_lower_bound() ) );

 mcf.chg_cost( -2 , 0 );
 check( mcf.get_valid_upper_bound( true ) == 42 , "conditional upper bound "
	"after chg_cost(): " + str( mcf.get_valid_upper_bound( true ) ) );
 check( mcf.get_valid_lower_bound() == -6 , "lower bound after chg_cost(): "
	+ str( mcf.get_valid_lower_bound() ) );

 mcf.chg_ucap( INF , 0 );
 check( mcf.get_valid_lower_bound() == -INF , "lower bound with a negative "
	"cost on an arc of infinite capacity" );

 // all capacities infinite, i.e., U empty
 MCFBlock mcf2;
 load( mcf2 , 0 , 0 , {} , { 1 , 0 , 1 , 0 , 1 } );
 check( mcf2.get_valid_lower_bound() == 0 , "no capacities: lower bound " +
	str( mcf2.get_valid_lower_bound() ) );
 check( mcf2.get_valid_upper_bound( true ) == INF , "no capacities: "
	"conditional upper bound " + str( mcf2.get_valid_upper_bound( true ) ) );

 // a deleted arc has no cost
 MCFBlock mcf3;
 load( mcf3 , 2 , 2 );
 mcf3.remove_arc( 3 );
 check( mcf3.get_valid_upper_bound( true ) == 30 , "deleted arc: "
	"conditional upper bound " + str( mcf3.get_valid_upper_bound( true ) ) );
 check( mcf3.get_valid_lower_bound() == 0 , "deleted arc: lower bound " +
	str( mcf3.get_valid_lower_bound() ) );
 }

/*--------------------------------------------------------------------------*/
/// random changes against a copy of the data, with a fixed seed

static void test_random( void )
{
 std::mt19937 rg( 20260926 );
 auto rnd = [ & ]( int lo , int hi ) {  // uniform in [ lo , hi ]
  return( std::uniform_int_distribution< int >( lo , hi )( rg ) );
  };
 auto rcap = [ & ]() {  // a capacity, infinite once in six
  return( rnd( 0 , 5 ) ? double( rnd( 0 , 9 ) ) : INF );
  };

 MCFBlock mcf;
 load( mcf , 2 , 6 , U0 , C0 , B0 , 1 );  // node 3 dynamic
 auto fs = attach( mcf );
 const Index NS = mcf.get_NStaticArcs();
 const Index MA = mcf.get_MaxNArcs();

 // the copy of the data
 Index m = 5;
 Subset sn( SN0 ) , en( EN0 );
 Vec_CNumber c( C0 );
 Vec_FNumber u( U0 ) , b( B0 );
 std::vector< bool > cls( 5 , false ) , del( 5 , false );
 sn.resize( MA ); en.resize( MA ); c.resize( MA ); u.resize( MA );
 cls.resize( MA ); del.resize( MA );

 auto names = [ & ]( Index nmax ) {  // a random subset of distinct names
  Subset nms;
  for( Index i = 0 ; i < nmax ; ++i )
   if( ! rnd( 0 , 2 ) )
    nms.push_back( i );
  std::shuffle( nms.begin() , nms.end() , rg );
  return( nms );
  };

 for( int step = 0 ; step < 400 ; ++step ) {
  const int op = rnd( 0 , 11 );
  std::string what = "random step " + std::to_string( step ) + " op " +
                     std::to_string( op );
  try {
   switch( op ) {
    case( 0 ): {  // costs of a Range
     Index f = rnd( 0 , m ) , s = rnd( f , m + 2 );
     Vec_CNumber nc( s - f );
     for( auto & v : nc ) v = rnd( -3 , 9 );
     mcf.chg_costs( std::span< const double >( nc ) , Range( f , s ) );
     for( Index a = f ; a < std::min( s , m ) ; ++a )
      if( ! del[ a ] ) c[ a ] = nc[ a - f ];
     break;
     }
    case( 1 ): {  // costs of a Subset
     auto nms = names( m );
     Vec_CNumber nc( nms.size() );
     for( auto & v : nc ) v = rnd( -3 , 9 );
     for( Index k = 0 ; k < nms.size() ; ++k )
      if( ! del[ nms[ k ] ] ) c[ nms[ k ] ] = nc[ k ];
     mcf.chg_costs( std::span< const double >( nc ) , std::move( nms ) );
     break;
     }
    case( 2 ): {  // capacities of a Range
     Index f = rnd( 0 , m ) , s = rnd( f , m + 2 );
     Vec_FNumber nu( s - f );
     for( auto & v : nu ) v = rcap();
     mcf.chg_ucaps( std::span< const double >( nu ) , Range( f , s ) );
     for( Index a = f ; a < std::min( s , m ) ; ++a )
      if( ! del[ a ] ) u[ a ] = nu[ a - f ];
     break;
     }
    case( 3 ): {  // capacities of a Subset
     auto nms = names( m );
     Vec_FNumber nu( nms.size() );
     for( auto & v : nu ) v = rcap();
     for( Index k = 0 ; k < nms.size() ; ++k )
      if( ! del[ nms[ k ] ] ) u[ nms[ k ] ] = nu[ k ];
     mcf.chg_ucaps( std::span< const double >( nu ) , std::move( nms ) );
     break;
     }
    case( 4 ): {  // deficits of a Range
     Index f = rnd( 0 , 4 ) , s = rnd( f , 5 );
     Vec_FNumber nb( s - f );
     for( auto & v : nb ) v = rnd( -5 , 5 );
     mcf.chg_dfcts( std::span< const double >( nb ) , Range( f , s ) );
     for( Index i = f ; i < std::min( s , Index( 4 ) ) ; ++i )
      b[ i ] = nb[ i - f ];
     break;
     }
    case( 5 ): {  // deficits of a Subset
     auto nms = names( 4 );
     Vec_FNumber nb( nms.size() );
     for( auto & v : nb ) v = rnd( -5 , 5 );
     for( Index k = 0 ; k < nms.size() ; ++k )
      b[ nms[ k ] ] = nb[ k ];
     mcf.chg_dfcts( std::span< const double >( nb ) , std::move( nms ) );
     break;
     }
    case( 6 ): {  // close a Range
     Index f = rnd( 0 , m ) , s = rnd( f , m + 2 );
     mcf.close_arcs( Range( f , s ) );
     for( Index a = f ; a < std::min( s , m ) ; ++a )
      if( ! del[ a ] ) cls[ a ] = true;
     break;
     }
    case( 7 ): {  // open a Subset
     auto nms = names( m );
     for( auto a : nms )
      if( ! del[ a ] ) cls[ a ] = false;
     mcf.open_arcs( std::move( nms ) );
     break;
     }
    case( 8 ): {  // close or open an arc
     if( ! m ) break;
     Index a = rnd( 0 , m - 1 );
     const bool close = rnd( 0 , 1 );
     if( close )
      mcf.close_arc( a );
     else
      mcf.open_arc( a );
     if( ! del[ a ] )
      cls[ a ] = close;
     break;
     }
    case( 9 ): {  // add an arc, no self-loop
     Index s = rnd( 1 , 4 ) , e = rnd( 1 , 3 );
     if( e >= s ) ++e;
     double cc = rnd( -3 , 9 ) , cu = rcap();
     Index a = NS;
     while( ( a < m ) && ( ! del[ a ] ) ) ++a;
     const Index got = mcf.add_arc( s , e , cc , cu );
     if( a >= MA ) {
      check( got == Inf< Index >() , what + ": add_arc() when full" );
      break;
      }
     check( got == a , what + ": add_arc() name " + std::to_string( got ) +
	    " instead of " + std::to_string( a ) );
     sn[ a ] = s; en[ a ] = e; c[ a ] = cc; u[ a ] = cu;
     cls[ a ] = del[ a ] = false;
     if( a == m ) ++m;
     break;
     }
    case( 10 ): {  // remove an arc
     if( m <= NS ) break;
     Index a = rnd( NS , m - 1 );
     mcf.remove_arc( a );
     if( del[ a ] ) break;
     cls[ a ] = false;
     if( a == m - 1 ) {
      --m;
      while( ( m > NS ) && del[ m - 1 ] ) { del[ m - 1 ] = false; --m; }
      }
     else
      del[ a ] = true;
     break;
     }
    default: {  // one cost, one capacity, one deficit
     if( m ) {
      Index a = rnd( 0 , m - 1 );
      double cc = rnd( -3 , 9 ) , cu = rcap();
      mcf.chg_cost( cc , a );
      mcf.chg_ucap( cu , a );
      if( ! del[ a ] ) c[ a ] = cc;
      u[ a ] = cu;
      }
     Index i = rnd( 0 , 3 );
     double nb = rnd( -5 , 5 );
     mcf.chg_dfct( nb , i );
     b[ i ] = nb;
     }
    }
   }
  catch( std::exception & ex ) {
   check( false , what + ": " + ex.what() );
   break;
   }

  auto md = take( fs );
  check_no_concern( md , what );

  check( mcf.get_NArcs() == m , what + ": " + std::to_string( mcf.get_NArcs() )
	 + " arcs instead of " + std::to_string( m ) );
  if( mcf.get_NArcs() != m )
   break;
  const int before = failures;
  for( Index a = 0 ; a < m ; ++a ) {
   check( mcf.is_deleted( a ) == del[ a ] , what + ": deletion of arc " +
	  std::to_string( a ) );
   if( del[ a ] )
    continue;
   check( ( mcf.get_SN( a ) == sn[ a ] ) && ( mcf.get_EN( a ) == en[ a ] ) ,
	  what + ": nodes of arc " + std::to_string( a ) );
   check( same( mcf.get_C( a ) , c[ a ] ) , what + ": cost of arc " +
	  std::to_string( a ) + " is " + str( mcf.get_C( a ) ) + " instead of "
	  + str( c[ a ] ) );
   check( mcf.get_U( a ) == u[ a ] , what + ": capacity of arc " +
	  std::to_string( a ) + " is " + str( mcf.get_U( a ) ) +
	  " instead of " + str( u[ a ] ) );
   check( mcf.is_closed( a ) == cls[ a ] , what + ": closure of arc " +
	  std::to_string( a ) );
   }
  for( Index i = 0 ; i < 4 ; ++i )
   check( mcf.get_B( i ) == b[ i ] , what + ": deficit of node " +
	  std::to_string( i ) );
  check_abstract( mcf , what );
  if( failures > before )  // the copy is no longer to be trusted
   break;
  }

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::vector< std::pair< std::string , void (*)( void ) > > tests = {
  { "load_no_data" , test_load_no_data } ,
  { "load_infinite_cost" , test_load_infinite_cost } ,
  { "load_errors" , test_load_errors } ,
  { "dmx_round_trip" , test_dmx_round_trip } ,
  { "netcdf_round_trip" , test_netcdf_round_trip } ,
  { "zero_arcs" , test_zero_arcs } ,
  { "costs" , test_costs } ,
  { "capacities" , test_capacities } ,
  { "deficits" , test_deficits } ,
  { "close_open" , test_close_open } ,
  { "mod_params" , test_mod_params } ,
  { "add_remove" , test_add_remove } ,
  { "remove_before_generate" , test_remove_before_generate } ,
  { "remove_all" , test_remove_all } ,
  { "range_ending_deleted" , test_range_ending_deleted } ,
  { "abstract_to_physical" , test_abstract_to_physical } ,
  { "flow_by_hand" , test_flow_by_hand } ,
  { "bounds" , test_bounds } ,
  { "random" , test_random } };

 for( auto & [ name , test ] : tests )
  try {
   test();
   }
  catch( std::exception & e ) {
   check( false , name + ": exception " + e.what() );
   }

 if( failures )
  std::cout << failures << " checks failed" << std::endl;
 else
  std::cout << "All tests passed!!" << std::endl;

 return( failures );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/

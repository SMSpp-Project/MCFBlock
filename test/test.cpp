/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing MCFBlock and MCFSolver.
 *
 * An instance of a MCF in DIMACS format is read from file in both an object
 * of class MCFC derived from MCFClass, and a MCFBlock to which a
 * MCFSolver<MCFC> is attached. The MCF problem is then repeatedly solved
 * with several changes in costs / capacities / deficits and arcs openings /
 * closures. The same operations are performed on the two solvers, and the
 * results are compared.
 *
 * \version 1.00
 *
 * \date 01 - 10 - 2018
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINES -----------------------------------*/
/*--------------------------------------------------------------------------*/
/* If any of the following macros is defined, then the corresponding
 * :MCFClass solver is included and the corresponding version of
 * MCFSolver<> can be tested.
 *
 * - HAVE_CSCL2      for the CS2 class
 *
 * - HAVE_CPLEX      for the MCFCplex class
 *
 * - HAVE_MFSMX      for the MCFSimplex class
 *
 * - HAVE_MFZIB      for the MCFZIB class
 *
 * - HAVE_RELAX      for the RelaxIV class
 *
 * - HAVE_SPTRE      for the MCFCplex class
 *
 * - HAVE_CPLEX      for the SPTree class; note that SPTree cannot solve
 *                   most MCF instances, except those with SPT structure
 *
 * Thus, the choice of the specific :MCFClass solver can be done in the
 * makefile with a simple -DHAVE_* argument to the compiler.
 */

#define NMS_IS_USED 0

// if NMS_IS_USED > 0, then the Chg****() routines are fed with a
// non-consecutive set of names; otherwise, all the involved arcs are
// consecutive

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef HAVE_CSCL2
 #include "CS2.h"
 #define MCFC CS2
#endif

#ifdef HAVE_CPLEX
 #include "MCFCplex.h"
 #define MCFC MCFCplex
#endif

#ifdef HAVE_MFSMX
 #include "MCFSimplex.h"
 #define MCFC MCFSimplex
#endif

#ifdef HAVE_MFZIB
 #include "MCFZIB.h"
 #define MCFC MCFZIB
#endif

#ifdef HAVE_RELAX
 #include "RelaxIV.h"
 #define MCFC RelaxIV
#endif

#ifdef HAVE_SPTRE
 #include "SPTree.h"
 #define MCFC SPTree
#endif

#include "MCFBlock.h"
#include "MCFSolver.h"
#include "FakeSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define solver_name( snm ) "MCFSolver<" #snm ">"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace MCFClass_di_unipi_it;
#else
 using namespace std;
#endif

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

MCFBlock * oMCFB;  // original MCFBlock
MCFBlock * dMCFB;  // "derived" (i.e., R3B) MCFBlock
MCFBlock * sMCFB;  // MCFBlock that is solved
MCFBlock * mMCFB;  // MCFBlock that is modified

MCFClass * mcf;    // the MCFClass object

vector<MCFClass::Index> open;  // names of open (and closed) arcs
int opened;                    // number of opened arcs

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/

static inline double rndfctr( void )
{
 // return a random number between 0.5 and 2, with 50% probability of being
 // < 1
 double fctr = drand48() - 0.5;
 return( fctr < 0 ? - fctr : fctr * 4 );
 }

/*--------------------------------------------------------------------------*/

static inline void CreateProb( int Optns )
{
 bool reoptmz = Optns & 1;
 Optns /= 2;

 mcf = nullptr;  // unknown solver, or the required solver is not
                 // available due to the macroes settings

 #ifdef HAVE_RELAX  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RelaxIV *rlx = new RelaxIV();
  #if( AUCTION )
   if( Optns )
    rlx->SetPar( RelaxIV::kAuction , MCFClass::kYes );
  #endif
  mcf = rlx;
  cout << "RelaxIV";
  assert( false );  // RelaxIV not fully supported yet
 #endif
 #ifdef HAVE_SPTRE  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mcf = new SPTree();
  cout << "SPTree";
  assert( false );  // SPTree not fully supported yet
 #endif
 #ifdef HAVE_CPLEX  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
  MCFCplex *cpx = new MCFCplex();
  if( Optns >= 0 )
   cpx->SetPar( CPX_PARAM_NETPPRIIND , Optns );
  mcf = cpx;
  cout << "MCFCplex";
  assert( false );  // MCFCplex not fully supported yet
 #endif
 #ifdef HAVE_MFZIB  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
  MCFZIB *zib = new MCFZIB();
  bool PrmlSmplx = Optns & 1;
  char Prcng;
  switch( Optns / 2 ) {
   case( 0 ): Prcng = char( MCFZIB::kDantzig ); break;
   case( 1 ): Prcng = char( MCFZIB::kFrstElA ); break;
   default:   Prcng = char( MCFZIB::kMltPrPr );
   }
  if( ( ! PrmlSmplx ) && ( Prcng == MCFZIB::kDantzig ) )
   Prcng = char( MCFZIB::kMltPrPr );
  zib->SetAlg( PrmlSmplx , Prcng );
  mcf = zib;
  cout << "MCFZIB";
  assert( false );  // MCFZIB not fully supported yet
 #endif
 #ifdef HAVE_CSCL2  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mcf = new CS2();
  cout << "CS2";
  assert( false );  // CS2 not fully supported yet
 #endif
 #ifdef HAVE_MFSMX  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
  MCFSimplex *mcfs = new MCFSimplex();
  bool PrmlSmplx = Optns & 1;
  char Prcng;
  switch( Optns / 2 ) {
   case( 0 ): Prcng = char( MCFSimplex::kDantzig ); break;
   case( 1 ): Prcng = char( MCFSimplex::kFirstEligibleArc ); break;
   default:   Prcng = char( MCFSimplex::kCandidateListPivot );
   }
  if( ( ! PrmlSmplx ) && ( Prcng == MCFSimplex::kDantzig ) )
   Prcng = char( MCFSimplex::kCandidateListPivot );
  mcf = mcfs;
  cout << "MCFSimplex";
 #endif

 if( ! reoptmz )
  mcf->SetPar( MCFClass::kReopt , MCFClass::kNo );

 if( ! mcf )
  throw( std::logic_error( "No MCFClass defined" ) );
 }

/*--------------------------------------------------------------------------*/

static inline void load( char * fn )
{
 ifstream iFile( fn );
 if( ! iFile ) {
  cerr << "Can't open input file " << fn << endl;
  exit( 1 );
  }

 try {
  mcf->LoadDMX( iFile );  // load the MCFClass

  iFile.clear();
  iFile.seekg( 0 );       // rewind the file

  iFile >> *oMCFB;        // load the MCFBlock

  }
 catch( exception &e ) {
  cerr << "MCFClass: " << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }

 iFile.close();

 // open[ 0 .. opened - 1 ] = names of open arcs
 // open[ opened ... m ] = names of closed arcs
 open.resize( opened = mcf->MCFm() );  // all arcs are open
 for( int i = 0 ; i < opened ; i++ )
  open[ i ] = i;

 }  // end( load )

/*--------------------------------------------------------------------------*/

static inline bool SolveMCF( void ) 
{
 try {
  // solve the MCFClass- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mcf->SolveMCF();
  cout << "MCFClass = ";
  switch( mcf->MCFGetStatus() ) {
   case( MCFClass::kOK ):         cout << mcf->MCFGetFO();
                                  break;
   case( MCFClass::kUnfeasible ): cout << "        +INF";
                                  break;
   case( MCFClass::kUnbounded ):  cout << "        -INF";
                                  break;
   default:                       cout << "      Error!";
   }

  cout << endl;

  // solve the MCFBlock- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // before actually solving, if modifications are not made on the MCFBlock
  // then exploit the FakeSolver to map them

  if( sMCFB != mMCFB ) {
   FakeSolver * fs = dynamic_cast<FakeSolver *>(
				  (mMCFB->get_registered_solvers()).front() );
   assert( ! fs );
   Lst_sp_Mod & modlist = fs->get_Modification_list();

   if( mMCFB == oMCFB ) {  // modify the original, solve the R3
    for( auto mod : modlist )
     oMCFB->map_forward_Modification( dMCFB , mod );
    }
   else {                  // modify the R3, solve the original
    for( auto mod : modlist )
     oMCFB->map_back_Modification( dMCFB , mod );
    }

   modlist.clear();  // clear the processed Modification
   }

  cout << "MCFBlock = ";
  Solver * slvr = (sMCFB->get_registered_solvers()).front();
  int rtrn = slvr->compute( false );
  if( ( rtrn > Solver::kUnEval ) && ( rtrn <= Solver::kOK ) ) {
   cout << slvr->get_ub() << endl;
   return( true );
   }

  switch( rtrn ) {
   case( Solver::kInfeasible ):   cout << "        +INF";
                                  break;
   case( Solver::kUnbounded ):    cout << "        -INF";
                                  break;
   default:                       cout << "      Error!";
   }

  cout << endl;

  return( false );
  }
 catch( exception &e ) {
  cerr << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 long int seed = 1;
 MCFClass::Index n_change = 10;
 MCFClass::Index n_repeat = 10;
 int optns = 0;
 int mode = 0;

 switch( argc ) {
  case( 7 ): Str2Sthg( argv[ 6 ] , seed );
  case( 6 ): Str2Sthg( argv[ 5 ] , n_change );
  case( 5 ): Str2Sthg( argv[ 4 ] , n_repeat );
  case( 4 ): Str2Sthg( argv[ 3 ] , optns );
  case( 3 ): Str2Sthg( argv[ 2 ] , mode );
  case( 2 ): break;
  default: cerr <<
	   "Usage: MCFSolve <input file> [mode optns #repeats #changes seed"
		<< endl <<
	   "       mode: 0 = only one, 1 = orig -> solve, 2 = solve -> orig"
		<< endl <<
	   "       optns: bit 0 = re-optimize, other bits MCF-specific" 
		<< endl <<
	   "              Relax   : > 0 uses Auction"
		<< endl <<
	   "              Cplex   : network pricing parameter"
		<< endl <<
	   "              ZIB     : 1st bit == 1 ==> primal +"
		<< endl <<
	   "                        0 = Dantzig, 2 = First Eligible, 4 = MPP"
	        << endl <<
	   "              Simplex : 1st bit == 1 ==> primal +"
		<< endl <<
	   "                        0 = Dantzig, 2 = First Eligible, 4 = MPP"
	        << endl;
	   return( 1 );
  }


 // construction and loading of the objects - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // construct MCFClass- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CreateProb( optns );

 //mcf->SetMCFTime();  // do timing

 // construct MCFBlocks/MCFSolvers- - - - - - - - - - - - - - - - - - - - - -

 oMCFB = dynamic_cast<MCFBlock *>( Block::new_Block( "MCFBlock" ) );
 assert( oMCFB );

 // load the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 load( argv[ 1 ] );

 if( mode ) {                     // also use an R3 MCFBlock = copy
  dMCFB = dynamic_cast<MCFBlock *>( oMCFB->get_R3_Block() );  // construct it
  assert( dMCFB );

  if( mode == 1 ) {               // modify the original, solve the R3
   mMCFB = oMCFB;
   sMCFB = dMCFB;
   }
  else {                          // modify the R3, solve the original
   mMCFB = dMCFB;
   sMCFB = oMCFB;
   }

  // attach a FakeSolver to the modified one to syphoon off Modification
  mMCFB->register_Solver( Solver::new_Solver( "FakeSolver" ) );
  }
 else {                    // just use one MCFBlock
  sMCFB = mMCFB = oMCFB;
  dMCFB = nullptr;
  }

 // attach a "true" MCFSolver to the one that is actually solved
 sMCFB->register_Solver( Solver::new_Solver( solver_name( MCFC ) ) );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFClass::cIndex n = mcf->MCFn();
 MCFClass::cIndex m = mcf->MCFm();

 cout << ", n = " << n << ", m = " << m << endl;
 if( n_change > m )
  n_change = m;

 // compute min/max cost & max deficit- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFClass::CNumber c_max = - OPTtypes_di_unipi_it::Inf<MCFClass::CNumber>();
                                                        // max cost
 MCFClass::CNumber c_min = - c_max;                     // min cost
 MCFClass::FNumber u_max = 0;                           // max capacity

 for( MCFClass::Index i = 0 ; i < m ; i++ ) {
  MCFClass::cCNumber ci = mcf->MCFCost( i );
  if( ci < c_min )
   c_min = ci;

  if( ci > c_max )
   c_max = ci;

  MCFClass::cFNumber ui = mcf->MCFUCap( i );
  if( ui > u_max )
   u_max = ui;
  }

 bool nzdfct = false;

 for( MCFClass::Index i = 0 ; i < n ; )
  if( mcf->MCFDfct( i++ ) > 0 ) {
   nzdfct = true;
   break;
   }

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << "First call:\t\t ";
 cout.setf( ios::scientific, ios::floatfield );
 cout << setprecision( 6 );

 SolveMCF();
 
 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now, for n_repeat times:
 // - up tp n_change costs are changed, then the two problems are re-solved;
 // - up to n_change capacities are changed, then the two problems are
 //   re-solved, then the original capacities are restored;
 // - if the problem is not a circulation problem, 2 deficits are modified
 //   (adding and subtracting the same number), then the two problems are
 //   re-solved, then the original deficits are restored;
 // - up to n_change arcs are closed, then the two problems are re-solved;
 //   the same arcs arcs are re-opened, then the two problems are re-solved

 srand48( seed );  // seed the pseudo-random number generator

 MCFBlock::Vec_CNumber newcsts( n_change );
 MCFBlock::Vec_FNumber newcaps( max( n_change , n ) );

 while( n_repeat-- ) {

  cout << "Changing: ";

  // change costs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( drand48() <= 0.2 ) {  // ... but only in 20% of the cases
   MCFBlock::Index tochange = min( double( 1 ) , drand48() * n_change );
   cout << tochange << " cost - ";

   if( tochange == 1 ) {
    MCFBlock::CNumber newcst = c_min +
                          MCFBlock::CNumber( drand48() * ( c_max - c_min ) );

    MCFBlock::Index arc = MCFBlock::Index( drand48() * ( m - 1 ) );

    mcf->ChgCost( arc , newcst );
    mMCFB->chg_cost( arc , newcst );
    }
   else {
    for( MCFBlock::Index i = 0 ; i < n_change ; i++ )
     newcsts[ i ] = c_min + MCFBlock::CNumber( drand48() * ( c_max - c_min ) );

    // in 50% of the cases do a ranged change, in the others a sparse change
    if( drand48() <= 0.5 ) {
     MCFBlock::Index strt = drand48() * ( m - tochange );
     MCFBlock::Index stp = strt + tochange;
     mcf->ChgCosts( newcsts.data() , nullptr , strt , stp );
     mMCFB->chg_costs( newcsts.begin() , strt , stp );
     }
    else {
     MCFBlock::Vec_Index nms( m + 1 );
     for( MCFBlock::Index i = 0 ; i < m ; i++ )
      nms[ i ] = i;

     for( MCFBlock::Index i = 0 ; i < tochange ; i++ )
      swap( nms[ i ] , nms[ i + drand48() * ( m - i ) ] );

     auto end = nms.begin() + tochange;
     sort( nms.begin() , end );
     *end = OPTtypes_di_unipi_it::Inf<MCFClass::Index>();
     mcf->ChgCosts( newcsts.data() , nms.data() );
     nms.resize( tochange );
     mMCFB->chg_costs( newcsts.begin() , std::move( nms ) , true );
     }
    }
   }  // end( if( change costs ) )

  // change capacities- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( drand48() <= 0.2 ) {  // ... but only in 20% of the cases
   MCFBlock::Index tochange = min( double( 1 ) , drand48() * n_change );
   cout << tochange << " capacity - ";

   if( tochange == 1 ) {
    MCFBlock::Index arc = MCFBlock::Index( drand48() * ( m - 1 ) );
    MCFBlock::CNumber newcap = mcf->MCFUCap( arc ) * rndfctr();
    mcf->ChgUCap( arc , newcap );
    mMCFB->chg_ucap( arc , newcap );
    }
   else
    // in 50% of the cases do a ranged change, in the others a sparse change
    if( drand48() <= 0.5 ) {
     MCFBlock::Index strt = drand48() * ( m - tochange );
     MCFBlock::Index stp = strt + tochange;
     for( MCFBlock::Index i = 0 ; i < n_change ; i++ )
      newcaps[ i ] = mcf->MCFUCap( i + strt ) * rndfctr();
     mcf->ChgUCaps( newcaps.data() , nullptr , strt , stp );
     mMCFB->chg_ucaps( newcaps.begin() , strt , stp );
     }
    else {
     MCFBlock::Vec_Index nms( m + 1 );
     for( MCFBlock::Index i = 0 ; i < m ; i++ )
      nms[ i ] = i;

     for( MCFBlock::Index i = 0 ; i < tochange ; i++ ) {
      swap( nms[ i ] , nms[ i + drand48() * ( m - i ) ] );
      newcaps[ i ] = mcf->MCFUCap( nms[ i ]  ) * rndfctr();
      }

     auto end = nms.begin() + tochange;
     sort( nms.begin() , end );
     *end = OPTtypes_di_unipi_it::Inf<MCFClass::Index>();
     mcf->ChgUCaps( newcsts.data() , nms.data() );
     nms.resize( tochange );
     mMCFB->chg_ucaps( newcsts.begin() , std::move( nms ) , true );
     }
   }  // end( if( change capacities ) )

  // change deficits- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( drand48() <= 0.2 ) {  // ... but only in 20% of the cases
   cout << " 2 deficits - ";

   MCFClass::Index posn;
   MCFClass::Index negn;
   MCFClass::FNumber posd;
   MCFClass::FNumber negd;

   if( nzdfct ) {  // if there are nonzero deficits
    mcf->MCFDfcts( newcaps.data() );

    do
     posn = MCFClass::Index( drand48() * n );  // select node with positive
    while( newcaps[ posn ] <= 0 );             // deficit (one must exist)
    posd = newcaps[ posn ];

    do
     negn = MCFClass::Index( drand48() * n );  // select node with negative
    while( newcaps[ negn ] >= 0 );             // deficit (one must exist)
    negd = newcaps[ negn ];
    }
   else {
    posn = MCFClass::Index( drand48() * n );   // just select at random
    negn = MCFClass::Index( drand48() * n );
    posd = negd = 0;
    }

   MCFClass::FNumber Dlt = u_max / 5;
   if( drand48() <= 0.5 ) {  // in 50% of cases up, in 50% of cases down
    posd += Dlt;
    negd -= Dlt;
    }
   else {
    Dlt = min( Dlt , max( max( posd , - negd ) / 2 , double( 1 ) ) );
    posd -= Dlt;
    negd += Dlt;
    }

   mcf->ChgDfct( posn , posd );
   mcf->ChgDfct( negn , negd );
   mMCFB->chg_dfct( posn , posd );
   mMCFB->chg_dfct( negn , negd );

   }  // end( change deficits )

  // closing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( drand48() <= 0.2 ) {  // ... but only in 20% of the cases
   MCFBlock::Index tochange = max( min( MCFBlock::Index( opened / 2 ) ,
					MCFBlock::Index( drand48() * n_change
							 ) ) ,
				   MCFBlock::Index( 1 ) );
                                   // and at most half of the open ones
   cout << tochange << " close - ";

   MCFBlock::Vec_Index nms( tochange );
   for( int i = 0 ; i < tochange ; i++ ) {
    MCFBlock::Index pos = drand48() * opened;
    MCFBlock::Index arc = open[ pos ];
    open[ pos ] = open[ --opened ];
    open[ opened ] = arc;
    nms[ i ] = arc;
    mcf->CloseArc( arc );
    }

   mMCFB->close_arcs( std::move( nms ) );
   }

  // re-opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( drand48() <= 0.2 ) {  // ... but only in 20% of the cases
   MCFBlock::Index tochange = min( MCFBlock::Index( ( m - opened ) / 2 ) ,
				   MCFBlock::Index( drand48() * n_change ) );
   // and at most half of the closed ones

   if( tochange ) {
    cout << tochange << " open - ";

    MCFBlock::Vec_Index nms( tochange );
    for( int i = 0 ; i < tochange ; i++ ) {
     MCFBlock::Index pos = drand48() * ( m - opened );
     MCFBlock::Index arc = open[ opened + pos ];
     open[ opened + pos ] = open[ opened ];
     open[ opened++ ] = arc;
     nms[ i ] = arc;
     mcf->OpenArc( arc );
     } 

    mMCFB->open_arcs( std::move( nms ) );
    }
   }

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // yet, if the problem is either unfeasible or unbounded, re-load it in
  // both MCFClass and MCFBlock

  if( ! SolveMCF() )
   load( argv[ 1 ] );

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 //double tu , ts;
 //mcf->TimeMCF( tu , ts );
 //cout << "Time: MCF = " << tu + ts << ", MCF2 = ";
 //mcf2->TimeMCF( tu , ts );
 //cout << tu + ts << endl;
 
 // destroy objects and vectors - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete dMCFB;
 delete oMCFB;
 delete mcf;

 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/

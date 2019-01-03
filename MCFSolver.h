/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MCFSolver class, implementing the Solver interface, in
 * particular in its CDASolver version, for Min-Cost Flow problems as set by
 * MCFBlock.
 *
 * This is only a relatively thin wrapper class around solvers under the
 * MCFClass interface. To avoid a pointer to an internal object, the class is
 * template over the underlying :MCFClass object, which implies that most of
 * the code is in the header file.
 *
 * \version 0.10
 *
 * \date 30 - 09 - 2018
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFSolver
 #define __MCFSolver  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"
#include "MCFBlock.h"
#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 #if( OPT_USE_NAMESPACES )
  using namespace MCFClass_di_unipi_it;
 #endif

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFSolver_CLASSES Classes in MCFSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// CDASolver for MCFBlock
/** The MCFSolver implements the Solver interface for Min-Cost Flow problems
 * described by a MCFBlock. Because the linear MCF problem is a Linear Program
 * it has a(n exact) dual, and therefore MCFSolver implements the CDASolver
 * interface for also giving out dual information.
 *
 * This is only a relatively thin wrapper class around solvers under the
 * MCFClass interface. To avoid a pointer to an internal object, the class is
 * template over the underlying :MCFClass object, which implies that most of
 * the code is in the header file. */

template< typename MCFC >
class MCFSolver : public CDASolver , private MCFC {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 /*
 kUnEval = 0     compute() has not been called yet

 kUnbounded = kUnEval + 1     the model is provably unbounded

 kInfeasible                  the model is provably infeasible

 kBothInfeasible = kInfeasible + 1     both primal and dual infeasible

 kOK = 7         successful compute()
                 Any return value between kUnEval (excluded) and kOK
		 (included) means that the object ran smoothly

 kStopTime = kOK + 1          stopped because of time limit

 kStopIter                    stopped because of iteration limit

 kError = 15     compute() stopped because of unrecoverable error
                 Any return value >= kError means that the object was
		  forced to stop due to some error, e.g. of numerical nature

 kLowPrecision = kError + 1   a solution found but not provably optimal
 */

/*--------------------------------------------------------------------------*/

 /*
 intMaxIter = 0     maximum iterations for the next call to solve()

 intMaxSol          maximum number of different solutions to report

 intLogVerb         "verbosity" of the log

 intMaxDSol         maximum number of different dual solutions

 intLastParCDAS     first allowed parameter value for derived classes
 */

/*--------------------------------------------------------------------------*/

 /*
 dblMaxTime = 0    maximum time for the next call to solve()

 dblRelAcc         relative accuracy for declaring a solution optimal

 dblAbsAcc          absolute accuracy for declaring a solution optimal

 dblUpCutOff        upper cutoff for stopping the algorithm

 dblLwCutOff        lower cutoff for stopping the algorithm

 dblRAccSol          maximum relative error in any reported solution

 dblAAccSol          maximum absolute error in any reported solution

 dblFAccSol          maximum constraint violation in any reported solution

 dblRAccDSol         maximum relative error in any dual solution

 dblAAccDSol         maximum absolute error in any dual solution

 dblFAccDSol         maximum absolute error in any dual solution

 dblLastParCDAS      first allowed parameter value for derived classes
 */

/*--------------------------------------------------------------------------*/

// typedef double OFValue;

/*@} -----------------------------------------------------------------------*/
/*----------------- CONSTRUCTING AND DESTRUCTING MCFSolver -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from MCFClass. */

 MCFSolver( void ) : CDASolver() , MCFC() {
  static_assert( std::is_base_of< MCFClass , MCFC >::value ,
                 "MCFSolver: MCFC must inherit from MCFClass" );
  }

/*--------------------------------------------------------------------------*/
 /// destructor: it has to release all the Modifications

 virtual ~MCFSolver() { }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, MCFSolver maps the parameters of [CDA]Solver
 *
 *  intMaxIter = 0    maximum iterations for the next call to solve()
 *  intMaxSol         maximum number of different solutions to report
 *  intLogVerb        "verbosity" of the log
 *  intMaxDSol        maximum number of different dual solutions
 *
 *  dblMaxTime = 0    maximum time for the next call to solve()
 *  dblRelAcc         relative accuracy for declaring a solution optimal
 *  dblAbsAcc         absolute accuracy for declaring a solution optimal
 *  dblUpCutOff       upper cutoff for stopping the algorithm
 *  dblLwCutOff       lower cutoff for stopping the algorithm
 *  dblRAccSol        maximum relative error in any reported solution
 *  dblAAccSol        maximum absolute error in any reported solution
 *  dblFAccSol        maximum constraint violation in any reported solution
 *  dblRAccDSol       maximum relative error in any dual solution
 *  dblAAccDSol       maximum absolute error in any dual solution
 *  dblFAccDSol       maximum absolute error in any dual solution
 *
 * into the parameter of MCFClass
 *
 * kMaxTime = 0       max time 
 * kMaxIter           max number of iteration
 * kEpsFlw            tolerance for flows
 * kEpsDfct           tolerance for deficits
 * kEpsCst            tolerance for costs
 *
 * It then "extends" them, using
 *
 *  intLastParCDAS    first allowed parameter value for derived classes
 *  dblLastParCDAS    first allowed parameter value for derived classes
 *
 * In particular, one now has
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * and any other paramater of specific :MCFClass following. This is done
 * via the two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl, with a negative entry meaning "there is no such
 * parameter in MCFSolver".
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve

 virtual void set_Block( Block * block ) override
 {
  if( block == f_Block )  // actually doing nothing
   return;                // cowardly and silently return

  Solver::set_Block( block );  // attach to the new Block

  if( block ) {  // this is not just resetting everything
   auto MCFB = dynamic_cast< MCFBlock * >( block );
   if( ! MCFB )
    throw( std::invalid_argument(
		         "MCFSolver:set_Block: block must be a MCFBlock" ) );

   // load the new MCFBlock into the :MCFClass object
   MCFC::LoadNet( MCFB->get_NNodes() , MCFB->get_NArcs() ,
		  MCFB->get_NNodes() , MCFB->get_NArcs() ,
		  MCFB->get_U().size() ? MCFB->get_U().data() : nullptr ,
		  MCFB->get_C().size() ? MCFB->get_C().data() : nullptr ,
		  MCFB->get_B().size() ? MCFB->get_B().data() : nullptr ,
		  MCFB->get_SN().data() , MCFB->get_EN().data() );
   MCFC::PreProcess();

   // TODO: maybe log it
   }
  }  // end( set_Block )

/*--------------------------------------------------------------------------*/
 // set the ostream for the Solver log
 // not really, MCFClass objects are remarkably silent
 //
 // virtual void set_log( std::ostream *log_stream = nullptr ) override;

/*--------------------------------------------------------------------------*/

 virtual void set_par( const idx_type par , const int value ) override
 {
  if( Solver_2_MCFClass_int[ par ] >= 0 )
   this->MCFC::SetPar( Solver_2_MCFClass_int[ par ] , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual void set_par( const idx_type par , const double value ) override
 {
  if( Solver_2_MCFClass_dbl[ par ] >= 0 )
   this->MCFC::SetPar( Solver_2_MCFClass_dbl[ par ] , value );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the MCF encoded in the MCFBlock

 virtual int compute( bool changedvars = true ) override
 {
  const static std::vector<int> MCFstatus_2_sol_type = {
   kUnEval , Solver::kOK , kStopTime , kInfeasible , Solver::kUnbounded ,
   Solver::kError };

  // first, process any outstanding Modification
  process_outstanding_Modification();

  // then (try to) solve the MCF
  this->MCFC::SolveMCF();

  // now give out the result: note that the vector MCFstatus_2_sol_type[]
  // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
  // hence the returned status has to be shifted by + 1
  return( MCFstatus_2_sol_type[ this->MCFC::MCFGetStatus() + 1 ] );
  }

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 virtual OFValue get_lb( void ) override { return( this->MCFC::MCFGetDFO() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual OFValue get_ub( void )  override { return( this->MCFC::MCFGetFO() ); }

/*--------------------------------------------------------------------------*/

 virtual bool has_var_solution( void ) override
 {
  switch( this->MCFC::MCFGetStatus() ) {
   case( MCFClass::kOK ):
   case( MCFClass::kUnbounded ): return( true );
   default:                      return( false );
   }
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool has_dual_solution( void ) override
 {
  switch( this->MCFGetStatus() ) {
   case( MCFClass::kOK ):
   case( MCFClass::kUnfeasible ): return( true );
   default:                       return( false );
   }
  }

/*--------------------------------------------------------------------------*/
/*
 virtual bool is_var_feasible( void ) override { return( true ); }

 virtual bool is_dual_feasible( void ) override { return( true ); }
*/
/*--------------------------------------------------------------------------*/

 virtual void get_var_solution( void ) override
 {
  auto MCFB = static_cast< MCFBlock * >( f_Block );
  if( ! MCFB )
   return;

  MCFBlock::Vec_FNumber X( MCFB->get_NArcs() );
  this->MCFGetX( X.data() );
  MCFB->set_x( X.begin() , X.end() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual void get_dual_solution( void ) override
 {
  auto MCFB = static_cast< MCFBlock * >( f_Block );
  if( ! MCFB )
   return;
  
  MCFBlock::Vec_CNumber Pi( MCFB->get_NNodes() );
  this->MCFGetPi( Pi.data() );
  MCFB->set_pi( Pi.begin() , Pi.end() );
  
  MCFBlock::Vec_FNumber RC( MCFB->get_NArcs() );
  this->MCFGetRC( RC.data() );
  MCFB->set_rc( RC.begin() , RC.end() );
  }

/*--------------------------------------------------------------------------*/

 virtual bool new_var_solution( void ) override
 {
  return( this->HaveNewX() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool new_dual_solution( void )  override
 {
  return( this->HaveNewPi() );
  }

/*--------------------------------------------------------------------------*/
/*
 virtual void set_unbounded_threshold( const OFValue thr ) override { }
*/

/*--------------------------------------------------------------------------*/

 virtual bool has_var_direction( void ) override { return( true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual bool has_dual_direction( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/

 virtual void get_var_direction( void ) override
 {
  throw( std::logic_error(
		    "MCFSolver::get_var_direction() not implemented yet" ) );

  // TODO: implement using MCFC::MCFGetUnbCycl()
  // anyway, unsure if any current :MCFClass properly implemente the latter
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual void get_dual_direction( void ) override
 {
  throw( std::logic_error(
		   "MCFSolver::get_dual_direction() not implemented yet" ) );

  // TODO: implement using MCFC::MCFGetUnfCut()
  // anyway, unsure if any current :MCFClass properly implemente the latter
  }

/*--------------------------------------------------------------------------*/
/*
 virtual bool new_var_direction( void ) override { return( false ); }

 virtual bool new_dual_direction( void ) override{ return( false ); }
*/
/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/

/*
 virtual bool is_dual_exact( void ) const override { return( true ); }
*/
 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the MCFSolver
 *
 * Each MCFSolver< MCFC > may have its own extra int / double parameters. If
 * this is the case, it will have to specialize the following methods to
 * handle them. The general definition just handles the case of the
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * extra (int) parameter and otherwise issues the method of the base
 * CDASolver class, which is OK for each MCFC that does *not* have any extra
 * parameter of the corresponding type (apart from that). The get_*_par()
 * methods exploit the same two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl as the set_*_par(), with a negative entry meaning
 * "there is no such parameter in MCFSolver".
 *  @{ */

 virtual idx_type get_num_int_par( void ) const override
 {
  return( CDASolver::get_num_int_par() + 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual idx_type get_num_dbl_par( void ) const override
 {
  return( CDASolver::get_num_dbl_par() );
  }

/*--------------------------------------------------------------------------*/
 
 virtual int get_dflt_int_par( const idx_type par ) const override
 {
  if( par == intLastParCDAS )
   return( MCFClass::kYes );
  else
   return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 virtual double get_dflt_dbl_par( const idx_type par ) const override
 {
  return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 virtual int get_int_par( const idx_type par ) const override
 {
  if( Solver_2_MCFClass_int[ par ] >= 0 ) {
   int val;
   this->GetPar( Solver_2_MCFClass_int[ par ] , val );
   return( val );
   }
  else
   return( get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 virtual double get_dbl_par( const idx_type par ) const override
 {
  if( Solver_2_MCFClass_dbl[ par ] >= 0 ) {
   double val;
   this->GetPar( Solver_2_MCFClass_dbl[ par ] , val );
   return( val );
   }
  else
   return( get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 virtual idx_type int_par_str2idx( const std::string & name ) const override
 {
  if( name == "kReopt" )
   return( intLastParCDAS );

  return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual idx_type dbl_par_str2idx( const std::string & name ) const override
 {
  return( CDASolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 virtual const std::string & int_par_idx2str( const idx_type idx )
  const override
 {
  static const std::string my_name = "kReopt";

  if( idx == intLastParCDAS )
   return( my_name );

  return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual const std::string & dbl_par_idx2str( const idx_type idx )
  const override
 {
  return( CDASolver::dbl_par_idx2str( idx ) );
  }

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /*
 virtual void add_Modification( sp_Mod &mod ) {
  v_mod.push_back( mod );
  }
 */

/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 const static std::vector<int> Solver_2_MCFClass_int;
 // the (static const) map between Solver int parameters and MCFClass ones

 const static std::vector<int> Solver_2_MCFClass_dbl;
 // the (static const) map between Solver int parameters and MCFClass ones

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class MCFSolver )

/*@}  end( group( Solver_CLASSES ) ) ---------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

template< class MCFC >
void MCFSolver< MCFC >::process_outstanding_Modification( void )
{
 // no-frills loop: do them in order, with no attempt at optimizing
 while( ! v_mod.empty() ) {
  auto mod = v_mod.front();  // pick (a reference to) the first Modification

  /* Use a Lambda to define a "guts" of the method that can be called
     recursively. Note the trick of defining the std::function object and
     "passing" it to the lambda, which allows recursive calls. Note the need
     to explicitly capture "this" to use fields/methods of the class. */

  auto MCFB = static_cast< MCFBlock * >( f_Block );

  std::function< void( sp_Mod )> guts_of_poM;
  guts_of_poM = [ this , & guts_of_poM , MCFB ]( sp_Mod mod ) {
   // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   /* This requires to patiently sift through the possible Modification types
      to find what this Modification exactly is, and call the appropriate
      method of MCFClass. */

   // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<GroupModification>( mod );
    if( tmod ) {
     for( const auto & submod : tmod->v_sub_Modifications )
      guts_of_poM( submod );

     return;
     }
    }

   // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   /* Note: in the following we can assume that C, B and U are nonempty. This
      is because they can be empty only if they are so when the object is
      loaded. But if a Modification has been issued they are no longer empty
      (a Modification changin nothing from the "empty" state is not issued).
      */
   {
    const auto tmod = std::dynamic_pointer_cast<MCFBlockRngdMod>( mod );
    if( tmod ) {
     switch( tmod->f_type ) {
      case( MCFBlockMod::eChgCost ):
       if( tmod->f_stop == tmod->f_strt + 1 )
        MCFC::ChgCost( tmod->f_strt , MCFB->get_C( tmod->f_strt ) );
       else
	MCFC::ChgCosts( MCFB->get_C().data() + tmod->f_strt , nullptr,
			tmod->f_strt , tmod->f_stop );
       break;

      case( MCFBlockMod::eChgCaps ):
       if( tmod->f_stop == tmod->f_strt + 1 )
        MCFC::ChgUCap( tmod->f_strt , MCFB->get_U( tmod->f_strt ) );
       else
	MCFC::ChgUCaps( MCFB->get_U().data() + tmod->f_strt , nullptr,
			tmod->f_strt , tmod->f_stop );
       break;

      case( MCFBlockMod::eChgDfct ):
       if( tmod->f_stop == tmod->f_strt + 1 )
        MCFC::ChgDfct( tmod->f_strt , MCFB->get_B( tmod->f_strt ) );
       else
	MCFC::ChgDfcts( MCFB->get_B().data() + tmod->f_strt , nullptr,
			tmod->f_strt , tmod->f_stop );
       break;

      case( MCFBlockMod::eOpenArc ):
       for( auto arc = tmod->f_strt ; arc < tmod->f_stop ; )
	MCFC::OpenArc( arc++ );
       break;

      case( MCFBlockMod::eCloseArc ):
       for( auto arc = tmod->f_strt ; arc < tmod->f_stop ; )
	MCFC::CloseArc( arc++ );
       break;

      default:
       throw( std::invalid_argument( "unknown MCFBlockRngdMod type" ) );
      }

     return;
     }
    }

   // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   {
    const auto tmod = std::dynamic_pointer_cast<MCFBlockSbstMod>( mod );
    if( tmod ) {
     switch( tmod->f_type ) {
      case( MCFBlockMod::eOpenArc ):
       for( auto arc : tmod->f_nms )
	MCFC::OpenArc( arc );
       return;

      case( MCFBlockMod::eCloseArc ):
       for( auto arc : tmod->f_nms )
	MCFC::CloseArc( arc );
       return;
      }

     // have to InINF-terminate the vector of indices (damn!)
     MCFBlock::Vec_Index nmsI( tmod->f_nms.size() + 1 );
     *copy( tmod->f_nms.begin() , tmod->f_nms.end() , nmsI.begin() ) =
                                                       Inf<MCFBlock::Index>();
     switch( tmod->f_type ) {
      case( MCFBlockMod::eChgCost ): {
       MCFBlock::Vec_CNumber NCost( tmod->f_nms.size() );
       auto C = MCFB->get_C();
       for( MCFBlock::Index i = 0 ; i < NCost.size() ; i++ )
	NCost[ i ] = C[ nmsI[ i ] ];

       MCFC::ChgCosts( NCost.data() , nmsI.data() );
       break;
       }

      case( MCFBlockMod::eChgCaps ): {
       MCFBlock::Vec_FNumber NCap( tmod->f_nms.size() );
       auto U = MCFB->get_U();
       for( MCFBlock::Index i = 0 ; i < NCap.size() ; i++ )
	NCap[ i ] = U[ nmsI[ i ] ];

       MCFC::ChgUCaps( NCap.data() , nmsI.data() );
       break;
       }

      case( MCFBlockMod::eChgDfct ): {
       MCFBlock::Vec_FNumber NDfct( tmod->f_nms.size() );
       auto B = MCFB->get_B();
       for( MCFBlock::Index i = 0 ; i < NDfct.size() ; i++ )
	NDfct[ i ] = B[ nmsI[ i ] ];

       MCFC::ChgDfcts( NDfct.data() , nmsI.data() );
       break;
       }

      default:
       throw( std::invalid_argument( "unknown MCFBlockSbstMod type" ) );
      }

     return;
     }
    }

   // MCFBlockMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // note: this is checked after the previous two because they derive from
   // MCFBlockMod, and hence the std::dynamic_pointer_cast<> would suceed
   /*
   {
    const auto tmod = std::dynamic_pointer_cast<MCFBlockMod>( mod );
    if( tmod ) {
      this is the "nuclear option": the MCFBlock has been re-loaded
     MCFC::LoadNet( MCFB->get_NNodes() , MCFB->get_NArcs() ,
		    MCFB->get_NNodes() , MCFB->get_NArcs() ,
		    MCFB->get_U().size() ? MCFB->get_U().data() : nullptr ,
		    MCFB->get_C().size() ? MCFB->get_C().data() : nullptr ,
		    MCFB->get_B().size() ? MCFB->get_B().data() : nullptr ,
		    MCFB->get_SN().data() , MCFB->get_EN().data() );
     MCFC::PreProcess();
     return;
     }
    }
   */
   {
    const auto tmod = std::dynamic_pointer_cast<BlockMod>( mod );
    if( tmod ) {
     if( tmod->f_type == BlockMod::eReSetAll ) {
      // this is the "nuclear option": the MCFBlock has been re-loaded
      MCFC::LoadNet( MCFB->get_NNodes() , MCFB->get_NArcs() ,
		     MCFB->get_NNodes() , MCFB->get_NArcs() ,
		     MCFB->get_U().size() ? MCFB->get_U().data() : nullptr ,
		     MCFB->get_C().size() ? MCFB->get_C().data() : nullptr ,
		     MCFB->get_B().size() ? MCFB->get_B().data() : nullptr ,
		     MCFB->get_SN().data() , MCFB->get_EN().data() );
      MCFC::PreProcess();
      return;
      }
     }
    }

   // IMPORTANT NOTE: any remaining Modification is plainly ignored. It must
   // be an "abstract" Modification, which this Solver does not need to look
   // at

   };  // end( guts_of_poM ) - - - - - - - - - - - - - - - - - - - - - - - - -
       //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // finally, call the "guts of" - - - - - - - - - - - - - - - - - - - - - - -

  guts_of_poM( mod );  // now the actual call

  v_mod.pop_front();   // now the Modification is processed: remove it
  
  }  // end( while( there are Modification ) )
 }  // end( MCFSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MCFSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/






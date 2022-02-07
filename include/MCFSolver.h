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
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni
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
 using namespace MCFClass_di_unipi_it;

 class MCFSolverState;  // forward declaration of MCFSolverState

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

/** @} ---------------------------------------------------------------------*/
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

/** @} ---------------------------------------------------------------------*/
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

 void set_Block( Block * block ) override
 {
  if( block == f_Block )  // actually doing nothing
   return;                // cowardly and silently return

  Solver::set_Block( block );  // attach to the new Block

  if( block ) {  // this is not just resetting everything
   auto MCFB = dynamic_cast< MCFBlock * >( block );
   if( ! MCFB )
    throw( std::invalid_argument(
		         "MCFSolver:set_Block: block must be a MCFBlock" ) );

   bool owned = MCFB->is_owned_by( f_id );
   if( ( ! owned ) && ( ! MCFB->read_lock() ) )
    throw( std::logic_error( "cannot acquire read_lock on MCFBlock" ) );

   // load the new MCFBlock into the :MCFClass object
   MCFC::LoadNet( MCFB->get_MaxNNodes() , MCFB->get_MaxNArcs() ,
		  MCFB->get_NNodes() , MCFB->get_NArcs() ,
		  MCFB->get_U().empty() ? nullptr : MCFB->get_U().data() ,
		  MCFB->get_C().empty() ? nullptr : MCFB->get_C().data() ,
		  MCFB->get_B().empty() ? nullptr : MCFB->get_B().data() ,
		  MCFB->get_SN().data() , MCFB->get_EN().data() );
   // TODO: PreProcess() changes the internal data of the MCFSolver using
   //       information about how the data of the MCF is *now*. If the data
   //       changes, some of the deductions (say, reducing the capacity but
   //       of some arcs) may no longer be correct and they should be undone,
   //       there isn't any proper way to handle this. Thus, PreProcess() has
   //       to be disabled for now; maybe later on someone will take care to
   //       make this work (or maybe not).
   // MCFC::PreProcess();

   // once done, read_unlock the MCFBlock (if it was read-lock()-ed)
   if( ! owned )
    MCFB->read_unlock();

   // TODO: maybe log it
   }
  }  // end( set_Block )

/*--------------------------------------------------------------------------*/
 // set the ostream for the Solver log
 // not really, MCFClass objects are remarkably silent
 //
 // virtual void set_log( std::ostream *log_stream = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void set_par( idx_type par , int value ) override {
  if( Solver_2_MCFClass_int[ par ] >= 0 )
   this->MCFC::SetPar( Solver_2_MCFClass_int[ par ] , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void set_par( idx_type par , double value ) override {
  if( Solver_2_MCFClass_dbl[ par ] >= 0 )
   this->MCFC::SetPar( Solver_2_MCFClass_dbl[ par ] , value );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the MCF encoded in the MCFBlock

 int compute( bool changedvars = true ) override
 {
  const static std::vector< int > MCFstatus_2_sol_type = {
   kUnEval , Solver::kOK , kStopTime , kInfeasible , Solver::kUnbounded ,
   Solver::kError };

  if( ! f_Block )           // there is no [MCFBlock] to solve
   return( kBlockLocked );  // return error 

  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure
  
  // while [read_]locked, process any outstanding Modification
  process_outstanding_Modification();

  if( ! owned )             // if the [MCF]Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it

  // ensure the timer exists (or reset it)
  this->MCFC::SetMCFTime();

  // then (try to) solve the MCF
  this->MCFC::SolveMCF();

  // now give out the result: note that the vector MCFstatus_2_sol_type[]
  // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
  // hence the returned status has to be shifted by + 1
  return( MCFstatus_2_sol_type[ this->MCFC::MCFGetStatus() + 1 ] );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 double get_elapsed_time( void ) const override {
  return( this->MCFC::TimeMCF() );
  }
 
/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override { return( this->MCFC::MCFGetDFO() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void )  override { return( this->MCFC::MCFGetFO() ); }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override {
  switch( this->MCFC::MCFGetStatus() ) {
   case( MCFClass::kOK ):
   case( MCFClass::kUnbounded ): return( true );
   default:                      return( false );
   }
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_dual_solution( void ) override {
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
 /// write the "current" flow in the x ColVariable of the MCFBlock
 /** Write the "current" flow in the x ColVariable of the MCFBlock. To keep
  * the same format as MCFBlock::get_Solution() and
  * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
  * be used to "partly" save it. In particular, if solc != nullptr, it is
  * a SimpleConfiguration<int>, and solc->f_value == 2, then *nothing is
  * done*, since the Configuration is meant to say "only save/map the dual
  * solution". In all other cases, the flow solution is saved. */

 void get_var_solution( Configuration * solc = nullptr ) override
 {
  if( ! f_Block )  // no [MCF]Block to write to
   return;         // cowardly and silently return

  auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  auto MCFB = static_cast< MCFBlock * >( f_Block );
  MCFBlock::Vec_FNumber X( MCFB->get_NArcs() );
  this->MCFGetX( X.data() );
  MCFB->set_x( X.begin() , X.end() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution in the Constraint of the MCFBlock
 /** Write the "current" dual solution, i.e., node potentials and flow
  * reduced costs, in the dual variables of the Constraint (respectively,
  * the flow conservation constraints and bound ones) of the MCFBlock. To
  * keep the same format as MCFBlock::get_Solution() and
  * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
  * be used to "partly" save it. In particular, if solc != nullptr, it is
  * a SimpleConfiguration<int>, and solc->f_value == 1, then *nothing is
  * done*, since the Configuration is meant to say "only save/map the primal
  * solution". In all other cases, the flow solution is saved. */

 void get_dual_solution( Configuration * solc = nullptr ) override
 {
  if( ! f_Block )  // no [MCF]Block to write to
   return;         // cowardly and silently return

  auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );
  if( tsolc && ( tsolc->f_value == 1 ) )
   return;

  auto MCFB = static_cast< MCFBlock * >( f_Block );
  MCFBlock::Vec_CNumber Pi( MCFB->get_NNodes() );
  this->MCFGetPi( Pi.data() );
  MCFB->set_pi( Pi.begin() , Pi.end() );
  
  MCFBlock::Vec_FNumber RC( MCFB->get_NArcs() );
  this->MCFGetRC( RC.data() );
  MCFB->set_rc( RC.begin() , RC.end() );
  }

/*--------------------------------------------------------------------------*/

 bool new_var_solution( void ) override { return( this->HaveNewX() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool new_dual_solution( void )  override { return( this->HaveNewPi() ); }

/*--------------------------------------------------------------------------*/
/*
 virtual void set_unbounded_threshold( const OFValue thr ) override { }
*/

/*--------------------------------------------------------------------------*/

 bool has_var_direction( void ) override { return( true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_dual_direction( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
 /// write the current direction in the x ColVariable of the MCFBlock
 /** Write the unbounded primal direction, i.e., augmenting cycle with
  * negative cost and unbounded capacity, in the x ColVariable of the
  * MCFBlock. To keep the same format as MCFBlock::get_Solution() and
  * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
  * be used to "partly" save it. In particular, if solc != nullptr, it is
  * a SimpleConfiguration<int>, and solc->f_value == 2, then *nothing is done*,
  * since the Configuration is meant to say "only save/map the dual
  * information". In all other cases, the direction (cycle) is saved.
  *
  * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

 void get_var_direction( Configuration * dirc = nullptr ) override
 {
  auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( dirc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  throw( std::logic_error(
		    "MCFSolver::get_var_direction() not implemented yet" ) );

  // TODO: implement using MCFC::MCFGetUnbCycl()
  // anyway, unsure if any current :MCFClass properly implemente the latter
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the current dual direction in the Constraint of the MCFBlock
 /** Write the current unbounded dual direction, i.e., a cut separating two
  * shores to that the residual demand in one is greater than the capacity
  * across them, in the Constraint of the Block, in particular in the dual
  * variables of the flow conservation ones. To keep the same format as
  * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
  * the Configuration *solc can be used to "partly" save it. In particular, if
  * solc != nullptr, it is a SimpleConfiguration<int>, and solc->f_value == 1,
  * then *nothing is done*, since the Configuration is meant to say "only
  * save/map the primal information". In all other cases, the direction (cut)
  * is saved.
  *
  * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

 void get_dual_direction( Configuration * dirc = nullptr ) override
 {
  auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( dirc );
  if( tsolc && ( tsolc->f_value == 1 ) )
   return;

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
/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/

/*
 virtual bool is_dual_exact( void ) const override { return( true ); }
*/

/*--------------------------------------------------------------------------*/
 /// "publicize" MCFClass::WriteMCF
 /** Make the method
  *
  *      void WriteMCF( ostream &oStrm , int frmt = 0 )
  *
  * of the base (private) MCFClass public, so that it can be freely used. */

 using MCFC::WriteMCF;

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

 idx_type get_num_int_par( void ) const override {
  return( CDASolver::get_num_int_par() + 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( CDASolver::get_num_dbl_par() );
  }

/*--------------------------------------------------------------------------*/
 
 int get_dflt_int_par( idx_type par ) const override {
  if( par == intLastParCDAS )
   return( MCFClass::kYes );

  return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dflt_dbl_par( idx_type par ) const override {
  return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_int_par( idx_type par ) const override {
  if( Solver_2_MCFClass_int[ par ] >= 0 ) {
   int val;
   this->GetPar( Solver_2_MCFClass_int[ par ] , val );
   return( val );
   }

  return( get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dbl_par( idx_type par ) const override {
  if( Solver_2_MCFClass_dbl[ par ] >= 0 ) {
   double val;
   this->GetPar( Solver_2_MCFClass_dbl[ par ] , val );
   return( val );
   }

  return( get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 idx_type int_par_str2idx( const std::string & name ) const override {
  if( name == "kReopt" )
   return( intLastParCDAS );

  return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  return( CDASolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 const std::string & int_par_idx2str( idx_type idx ) const override {
  static const std::string my_name = "kReopt";

  if( idx == intLastParCDAS )
   return( my_name );

  return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( idx_type idx ) const override {
  return( CDASolver::dbl_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR HANDLING THE State OF THE MCFSolver -------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the State of the MCFSolver
 *  @{ */

 State * get_State( void ) const override;

/*--------------------------------------------------------------------------*/

 void put_State( const State & state ) override;

/*--------------------------------------------------------------------------*/

 void put_State( State && state ) override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /** The only reason why MCFSolver::add_Modification() needs be defined is to
  * properly react to NBModification. Indeed, the correct reaction is to
  * *immediately* reload the MCFBlock, besides clearing the list of
  * Modification as Solver::add_Modification() already does. The issue is
  * that if arcs/nodes are added/deleted after the NBModification is issued
  * but before it is processed, then the number of nodes/arcs at the moment
  * in which the NBModification is processed is different from that at the
  * moment in which is is issued, which may break the "naming convention"
  * (because the name of, say, a newly created arc depends on the current
  * state and/or number of the arcs).
  *
  * Important note: THIS VERSION ONLY WORKS PROPERLY IF THE MCFBlock IS
  * "FRESHLY MINTED", I.E., THERE ARE NO CLOSED OR DELETED ARCS.
  *
  * This should ordinarily always happen, as whenever the MCFBlock is changed
  * the NBModification is immediately issued. The problem may come if the
  * MCFBlock is a R3Block of another MCFBlock which is loaded and then
  * further modified, and the NBModification to this MCFBlock is generated by
  * a map_forward_Modification() of the NBModification to the original
  * MCFBlock: then, this MCFBlock may be copied from a MCFBlock that has
  * closed or deleted arcs and this method would not work. */

 void add_Modification( sp_Mod &mod ) override
 {
  if( const auto tmod = std::dynamic_pointer_cast<NBModification>( mod ) ) {
   // this is the "nuclear option": the MCFBlock has been re-loaded, so
   // the MCFClass solver also has to (immediately)
   auto MCFB = static_cast< MCFBlock * >( f_Block );
   MCFC::LoadNet( MCFB->get_MaxNNodes() , MCFB->get_MaxNArcs() ,
		  MCFB->get_NNodes() , MCFB->get_NArcs() ,
		  MCFB->get_U().empty() ? nullptr : MCFB->get_U().data() ,
		  MCFB->get_C().empty() ? nullptr : MCFB->get_C().data() ,
		  MCFB->get_B().empty() ? nullptr : MCFB->get_B().data() ,
		  MCFB->get_SN().data() , MCFB->get_EN().data() );
   // TODO: PreProcess() changes the internal data of the MCFSolver using
   //       information about how the data of the MCF is *now*. If the
   //       data changes, some of the deductions (say, reducing the capacity
   //       of some arcs) may no longer be correct and they should be undone,
   //       but there isn't any proper way to handle this. Thus, PreProcess()
   //       has to be disabled for now; maybe later on someone will take
   //       care to make this work (or maybe not).
   // MCFC::PreProcess();
   // besides, any outstanding modification makes no sense any longer
   v_mod.clear();
   }
  else
   v_mod.push_back( mod );
  }

/*--------------------------------------------------------------------------*/
/*-------------------------------- FRIENDS ---------------------------------*/
/*--------------------------------------------------------------------------*/

 friend class MCFSolverState;  // make MCFSolverState friend

/** @} ---------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MCFSolverState ---------------------------*/
/*--------------------------------------------------------------------------*/
/// class to describe the "internal state" of a MCFSolver
/** Derived class from State to describe the "internal state" of a MCFSolver,
 *  i.e., a MCFClass::MCFState (*). Since MCFClass::MCFState does not allow
 *  serialization, all that part does not work.  */

class MCFSolverState : public State {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------- CONSTRUCTING AND DESTRUCTING MCFSolverState ----------------*/

 /// constructor, doing everything or nothing.
 /** Constructor of MCFSolverState. If provided with a pointer to a MCFSolver
  * it immediately copies its "internal state", which is the only way in which
  * the MCFSolverState can be initialised out of an existing MCFSolver. If
  * nullptr is passed (as by default), then an "empty" MCFSolverState is
  * constructed that can only be filled by calling deserialize().
  *
  * Note: to avoid having to duplicate the SMSpp_insert_in_factory_cpp call
  *       for every MCFClass, the pointer is directly that of a MCFClass,
  *       since every MCFSolver derives from a :MCFClass and we only need
  *       access to MCFGetState(). */

 MCFSolverState( MCFClass * mcfc = nullptr ) : State() {
  f_state = mcfc ? mcfc->MCFGetState() : nullptr;
  }

/*--------------------------------------------------------------------------*/
 /// de-serialize a MCFSolverState out of netCDF::NcGroup
 /** Should de-serialize a MCFSolverState out of netCDF::NcGroup, but in
  * fact it does not work. */

 void deserialize( const netCDF::NcGroup & group ) override {
  f_state = nullptr;
  }

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~MCFSolverState() { delete f_state; }

/*---------- METHODS DESCRIBING THE BEHAVIOR OF A MCFSolverState -----------*/

 /// serialize a MCFSolverState into a netCDF::NcGroup
 /** The method should serialize the MCFSolverState into the provided
  * netCDF::NcGroup, so that it can later be read back by deserialize(), but
  * in fact it does not work.*/

 void serialize( netCDF::NcGroup & group ) const override {}

/*-------------------------------- FRIENDS ---------------------------------*/

 template< class MCFC >
 friend class MCFSolver;  // make MCFSolver friend

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "MCFSolverState [" << this << "]";
  }

/*--------------------------- PROTECTED FIELDS -----------------------------*/

 MCFClass::MCFStatePtr f_state;   ///< the (pointer to) MCFState

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( MCFSolverState ) )

/** @} end( group( MCFSolver_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

template< class MCFC >
State * MCFSolver< MCFC >::get_State( void ) const {
 return( new MCFSolverState( const_cast< MCFSolver< MCFC > * >( this ) ) );
 }

/*--------------------------------------------------------------------------*/

template< class MCFC >
void MCFSolver< MCFC >::put_State( const State & state ) {
 // if state is not a const MCFSolverState &, exception will be thrown
 auto s = dynamic_cast< const MCFSolverState & >( state );

 this->MCFPutState( s.f_state );
 }

/*--------------------------------------------------------------------------*/

template< class MCFC >
void MCFSolver< MCFC >::put_State( State && state ) {
 // if state is not a MCFSolverState &&, exception will be thrown
 auto s = dynamic_cast< MCFSolverState && >( state );

 this->MCFPutState( s.f_state );
 }

/*--------------------------------------------------------------------------*/

template< class MCFC >
void MCFSolver< MCFC >::process_outstanding_Modification( void )
{
 // no-frills loop: do them in order, with no attempt at optimizing

 for( auto mod = front() ; mod ; mod = front() ) {
  /* Use a Lambda to define a "guts" of the method that can be called
     recursively. Note the trick of defining the std::function object and
     "passing" it to the lambda, which allows recursive calls. Note the need
     to explicitly capture "this" to use fields/methods of the class. */

  auto MCFB = static_cast< MCFBlock * >( f_Block );

  std::function< void( c_p_Mod )> guts_of_poM;
  guts_of_poM = [ this , & guts_of_poM , MCFB ]( c_p_Mod mod ) {
   // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   /* This requires to patiently sift through the possible Modification types
      to find what this Modification exactly is, and call the appropriate
      method of MCFClass. */

   // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( const auto tmod = dynamic_cast< GroupModification * const >( mod ) ) {
    for( const auto & submod : tmod->sub_Modifications() )
     guts_of_poM( submod.get() );

    return;
    }


   // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   /* Note: in the following we can assume that C, B and U are nonempty. This
      is because they can be empty only if they are so when the object is
      loaded. But if a Modification has been issued they are no longer empty
      (a Modification changin nothing from the "empty" state is not issued).
      */
   if( const auto tmod = dynamic_cast< MCFBlockRngdMod * const >( mod ) ) {
    switch( tmod->type() ) {
     case( MCFBlockMod::eChgCost ):
      if( tmod->rng().second == tmod->rng().first + 1 )
       MCFC::ChgCost( tmod->rng().first , MCFB->get_C( tmod->rng().first ) );
      else
       MCFC::ChgCosts( MCFB->get_C().data() + tmod->rng().first , nullptr,
		       tmod->rng().first , tmod->rng().second );
      break;

     case( MCFBlockMod::eChgCaps ):
      if( tmod->rng().second == tmod->rng().first + 1 )
       MCFC::ChgUCap( tmod->rng().first , MCFB->get_U( tmod->rng().first ) );
      else
       MCFC::ChgUCaps( MCFB->get_U().data() + tmod->rng().first , nullptr,
		       tmod->rng().first , tmod->rng().second );
      break;

     case( MCFBlockMod::eChgDfct ):
      if( tmod->rng().second == tmod->rng().first + 1 )
       MCFC::ChgDfct( tmod->rng().first , MCFB->get_B( tmod->rng().first ) );
      else
       MCFC::ChgDfcts( MCFB->get_B().data() + tmod->rng().first , nullptr,
		       tmod->rng().first , tmod->rng().second );
      break;

     case( MCFBlockMod::eOpenArc ):
      for( auto arc = tmod->rng().first ; arc < tmod->rng().second ; )
       MCFC::OpenArc( arc++ );
      break;

     case( MCFBlockMod::eCloseArc ):
      for( auto arc = tmod->rng().first ; arc < tmod->rng().second ; )
       MCFC::CloseArc( arc++ );
      break;

     case( MCFBlockMod::eAddArc ): {
      auto arc = MCFC::AddArc( MCFB->get_SN( tmod->rng().first ) ,
			       MCFB->get_EN( tmod->rng().first ) ,
			       MCFB->get_U( tmod->rng().first ) ,
			       MCFB->get_C( tmod->rng().first ) );
      if( arc != tmod->rng().first )
       throw( std::logic_error( "name mismatch in AddArc()" ) );
      break;
      }
     case( MCFBlockMod::eRmvArc ):
      MCFC::DelArc( tmod->rng().second - 1 );
      break;

     default:
      throw( std::invalid_argument( "unknown MCFBlockRngdMod type" ) );
     }

    return;
    }

   // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( const auto tmod = dynamic_cast< MCFBlockSbstMod * const >( mod ) ) {
    switch( tmod->type() ) {
     case( MCFBlockMod::eOpenArc ):
      for( auto arc : tmod->nms() )
       MCFC::OpenArc( arc );
      return;

     case( MCFBlockMod::eCloseArc ):
      for( auto arc : tmod->nms() )
       MCFC::CloseArc( arc );
      return;
     }

    // have to InINF-terminate the vector of indices (damn!)
    MCFBlock::Subset nmsI( tmod->nms().size() + 1 );
    *copy( tmod->nms().begin() , tmod->nms().end() , nmsI.begin() ) =
                                                       Inf<MCFBlock::Index>();
    switch( tmod->type() ) {
     case( MCFBlockMod::eChgCost ): {
      MCFBlock::Vec_CNumber NCost( tmod->nms().size() );
      auto C = MCFB->get_C();
      for( MCFBlock::Index i = 0 ; i < NCost.size() ; i++ )
       NCost[ i ] = C[ nmsI[ i ] ];

      MCFC::ChgCosts( NCost.data() , nmsI.data() );
      break;
      }

     case( MCFBlockMod::eChgCaps ): {
      MCFBlock::Vec_FNumber NCap( tmod->nms().size() );
      auto U = MCFB->get_U();
      for( MCFBlock::Index i = 0 ; i < NCap.size() ; i++ )
       NCap[ i ] = U[ nmsI[ i ] ];

      MCFC::ChgUCaps( NCap.data() , nmsI.data() );
      break;
      }

     case( MCFBlockMod::eChgDfct ): {
      MCFBlock::Vec_FNumber NDfct( tmod->nms().size() );
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

   // any remaining Modification is plainly ignored, since it must be an
   // "abstract" Modification, which this Solver does not need to look at

   };  // end( guts_of_poM ) - - - - - - - - - - - - - - - - - - - - - - - - -
       //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // finally, call the "guts of" - - - - - - - - - - - - - - - - - - - - - - -

  guts_of_poM( mod.get() );  // now the actual call

  pop_front();   // now the Modification is processed: remove it
  
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






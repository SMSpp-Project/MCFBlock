/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MCFBlock, which implements the Block
 * concept [see Block.h] for (linear) Min-Cost Flow problems.
 *
 * \version 0.11
 *
 * \date 23 - 02 - 2019
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

#ifndef __MCFBlock
 #define __MCFBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "Solution.h"
#include "LinearFunction.h"
#include "FRowConstraint.h"
#include "FRealObjective.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class MCFBlock;     // forward declaration of MCFBlock

 class MCFSolution;  // forward declaration of MCFSolution

/*--------------------------------------------------------------------------*/
/*----------------------- MCFBlock-RELATED TYPES ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFBlock_TYPES MCFBlock-related types
 *  @{ */

 typedef MCFBlock * p_MCFBlock;
 ///< a pointer to MCFBlock

 typedef std::vector<p_MCFBlock> Vec_MCFBlock;
 ///< a vector of pointers to MCFBlock

 typedef Vec_MCFBlock::iterator Vec_MCFBlock_it;
 ///< iterator for a Vec_MCFBlock

 typedef const std::vector<p_MCFBlock> c_Vec_MCFBlock;
 ///< a const vector of pointers to MCFBlock

 typedef c_Vec_MCFBlock::iterator c_Vec_MCFBlock_it;
 ///< iterator for a c_Vec_MCFBlock

/** @}  end( group( MCFBlock_TYPES ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFBlock_CLASSES Classes in MCFBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFBlock --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the (linear) Min-Cost Flow problem
/** The MCFBlock class implements the Block concept [see Block.h] for the
 * (linear) Min-Cost Flow (MCF) problem.
 *
 * The data of the problem consist of a (directed) graph G = ( N , A ) with
 * n = |N| nodes and m = |A| (directed) arcs. Each node `i' has a deficit
 * b[ i ], i.e., the amount of flow that is produced/consumed by the node:
 * source nodes (which produce flow) have negative deficits and sink nodes
 * (which consume flow) have positive deficits. Each arc `(i, j)' has an
 * upper capacity U[ i , j ] and a linear cost coefficient C[ i , j ]. Flow
 * variables X[ i , j ] represents the amount of flow to be sent on arc
 * (i, j). Parallel arcs, i.e., multiple copies of the same arc `(i, j)' are
 * in general allowed; it is expected that they have different costs (for
 * otherwise they can be merged into a unique arc), but this is not strictly
 * enforced. Multiple copies of some arc (i, j) can be seen as "total" flow
 * cost on that arc being a piecewise-linear convex function.
 *
 * The formulation of the problem is:
 * \f[
 *  \min \sum_{ (i, j) \in A } C[ i , j ] X[ i, j ]
 * \f]
 * \f[
 *  \sum_{ (j, i) \in A } X[ j , i ] -
 *  \sum_{ (i, j) \in A } X[ i , j ] = b[ i ] \quad i \in N      (1)
 * \f]
 * \f[
 *   0 \leq X[ i , j ] \leq U[ i , j ]   \quad (i, j) \in A      (2)
 * \f]
 * The n equations (1) are the flow conservation constraints and the 2m
 * inequalities (2) are the flow nonnegativity and capacity constraints.
 * At least one of the flow conservation constraints is redundant, as the
 * demands must be balanced (\f$\sum_{ i \in N } b[ i ] = 0\f$); indeed,
 * exactly n - ConnectedComponents( G ) flow conservation constraints are
 * redundant, as demands must be balanced in each connected component of G.
 *
 * The dual of the problem is:
 * \f[
 *  \max \sum_{ i \in N } Pi[ i ] b[ i ] -
 *       \sum_{ (i, j) \in A } W[ i , j ] U[ i , j ] -
 * \f]
 * \f[
 *  Pi[ j ] - Pi[ i ] - W[ i , j ] + Z[ i , j ] = C[ i , j ]
 *  \quad (i, j) \in A                                           (3)
 * \f]
 * \f[
 *  W[ i , j ] \geq 0 \quad (i, j) \in A                         (4)
 * \f]
 * \f[
 *  Z[ i , j ] \geq 0 \quad (i, j) \in A                         (5)
 * \f]
 *
 * Pi[] is said the vector of node potentials for the problem, W[] are bound
 * variables and Z[] are slack variables. Given Pi[], the quantities
 * \f[
 *  RC[ i , j ] =  C[ i , j ] + Q[ i , j ] * X[ i , j ] - Pi[ j ] + Pi[ i ]
 * \f]
 * are said the "reduced costs" of arcs, and are basically the dual variables
 * of the box constraints.
 *
 * A primal and dual feasible solution pair is optimal if and only if the
 * complementary slackness conditions
 * \f[
 *  RC[ i , j ] > 0 \Rightarrow X[ i , j ] = 0                   (6)
 * \f]
 * \f[
 *  RC[ i , j ] < 0 \Rightarrow X[ i , j ] = U[ i , j ]          (7)
 * \f]
 * are satisfied for all arcs (i, j) of A.
 */

class MCFBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * MCFBlock defines four main public types:
 *
 * - Index, the type of node indices;
 *
 * - Index, the type of arc indices;
 *
 * - FNumber, the type of flow variables, arc capacities, and node deficits;
 *
 * - CNumber, the type of flow costs, node potentials, and arc reduced costs;
 *
 * - FONumber, the type of objective function value.
 *
 * By re-defining the types in this section, some (but not all) solution
 * algorithms may be able to work with the "smallest" choice of data type 
 * that is capable of properly representing the data of the instances to be
 * solved. This may be relevant due to an important property of MCF problems:
 * *if all arc capacities and node deficits are integer, then there exists an
 * integral optimal primal solution*, and *if all arc costs are integer,
 * then there exists an integral optimal dual solution*. Even more
 * importantly, *many solution algorithms will in fact produce an integral
 * primal/dual solution for free*, because *every primal/dual solution they
 * generate during the solution process is naturally integral*. Therefore,
 * one can use integer data types to represent everything connected with
 * flows and/or costs if the corresponding data is integer in all instances
 * one needs to solve. This directly translates in significant memory savings
 * and/or speed improvements.
 *
 * However, while using a MCFBlock as a part of some larger problem, it may
 * be difficult to fully exploit this property: even if some Solver can
 * exploit it, not all ofthem may be able to (one example are Interior-Point
 * approaches, which require both flow and cost variables to be continuous),
 * and maybe some other aspects of the overall solution algorithm will require
 * general double data anyway. One should actually have Block template over
 * all these types to be able to fully exploit this property, which may be a
 * future evolution but is not what this implementation does. The current
 * choice is to use the "worst case scenario" where FNumber == CNumber ==
 * OFNumber == double, although the data types are left there and it is
 * therefore in principle possible to change this. Note, however, that the
 * above integrality property only holds for *linear* MCF problems. Should
 * the class be extended, by even allowing arc costs to be convex quadratic
 * (the simplest possible nonlionear extension), then a single arc with a
 * nonzero quadratic cost coefficient implies that optimal flows and
 * potentials may be fractional even if all the data of the problem
 * (comprised quadratic cost coefficients) is integer. Hence, for such a
 * setting FNumber == CNumber == OFNumber == double is actually *mandatory*,
 * for any reasonable algorithm will typically misbehave otherwise.
 @{ */

/*--------------------------------------------------------------------------*/

 typedef unsigned int Index;                 ///< index of a node / arc
 typedef const Index c_Index;                ///< a read-only Index

 typedef std::vector<Index> Vec_Index;       ///< a vector of Index
 typedef const Vec_Index c_Vec_Index;        ///< a const vector of Index

 typedef Vec_Index::iterator Vec_Index_it;   ///< iterator in Vec_Index
 typedef Vec_Index::const_iterator c_Vec_Index_it;
                                             ///< const iterator in Vec_Index

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double FNumber;                     ///< type of arc flow / deficit
 typedef const FNumber c_FNumber;            ///< a read-only FNumber

 typedef std::vector<FNumber> Vec_FNumber;   ///< a vector of FNumber
 typedef const Vec_FNumber c_Vec_FNumber;    ///< a const vector of FNumber

 typedef Vec_FNumber::iterator Vec_FNumber_it;   ///< iterator in Vec_FNumber
 typedef Vec_FNumber::const_iterator c_Vec_FNumber_it;
                                           ///< const iterator in Vec_FNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double CNumber;                     ///< type of arc cost / potential
 typedef const CNumber c_CNumber;            ///< a read-only CNumber

 typedef std::vector<CNumber> Vec_CNumber;   ///< a vector of CNumber
 typedef const Vec_CNumber c_Vec_CNumber;    ///< a const vector of CNumber

 typedef Vec_CNumber::iterator Vec_CNumber_it;   ///< iterator in Vec_CNumber
 typedef Vec_CNumber::const_iterator c_Vec_CNumber_it;
                                           ///< const iterator in Vec_CNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double FONumber; 
 /**< type of the objective function: has to hold sums of products of
    FNumber(s) by CNumber(s) */

 typedef const FONumber c_FONumber;            ///< a read-only FONumber

 typedef std::vector<FONumber> Vec_FONumber;   ///< a vector of FONumber
 typedef const Vec_FONumber c_Vec_FONumber;    ///< a const vector of FONumber

/*@} -----------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend MCFSolution;  ///< make MCFSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of MCFBlock, taking a pointer to the father (generic) Block
 /** Constructor of MCFBlock. It accepts a pointer to the father Block, which
  * can be of any type, defaulting to nullpt so that this can also be used as
  * the void constructor. */

 MCFBlock( Block *father = nullptr ) : Block( father ) , NNodes( 0 ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of MCFBlock: deletes the abstract representation, if any
 virtual ~MCFBlock() { guts_of_destructor(); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the MCF instance from memory
 /** Loads the MCF instance from memory. The parameters are what you expect:
  *
  * - n    is the number of nodes of the network
  *
  * - pSn  is the vector of the arc starting nodes;
  * - pEn  is the vector of the arc ending nodes;
  *
  * - pU   is the vector of the arc upper capacities; capacities must be
  *        nonnegative, but can be infinite; if pU is empty, then all
  *        capacities are taken to be infinite;
  *
  * - pC   is the vector of the arc costs; if pC is empty, then all arc costs
  *        are taken to be 0;
  *
  * - pB   is the vector of the node deficits; source nodes have negative
  *        deficits and sink nodes have positive deficits; if pB is empty,
  *        then all deficits are taken to be 0 (a circulation problem).
  *
  * The number m of arcs of the graph need not be explicitly provided because
  * it is the length of pEn and sEn; while pU and pC can be empty, they cannot
  * (unless the graph is empty of arcs). Conversely, n must be explicitly
  * provided because the only n-vector is pB, which can be empty.
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * MCFBlock then a NBModification (the "nuclear option") is issued. */

 virtual void load( c_Index n , c_Vec_Index & pEn , c_Vec_Index & pSn ,
		    c_Vec_FNumber & pU = {} , c_Vec_CNumber & pC = {} ,
		    c_Vec_FNumber & pB = {} );

/*--------------------------------------------------------------------------*/

 virtual void deserialize( netCDF::NcGroup&& group ,
			   Block *father = nullptr ) override;

/*--------------------------------------------------------------------------*/

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
  override final;

 /// generate the abstract variables of the MCF
 /** Method that generates the abstract variables of the MCF. These are the a
  * std::vector< ColVariable > with exactly m entries, the entry a = 0, ...,
  *  m - 1 corresponding to the flow on arc ( SN[ a ] , EN[ a ] ). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the MCF
 /** Method that generates the static constraint of the MCF. These are the:
  *
  * - the flow conservation equations, a std::vector<FRowConstraint> with
  *   exactly n entries, the entry i = 0, ..., n - 1 being the flow
  *   conservation of the node i;
  *
  * - the bound constraints, a std::vector< some derived class from
  *   OneVarConstraint > with exactly m entries, the entry a = 0, ..., m - 1
  *   being the bound constraints of the ColVariable x[ a ] corresponding to
  *   the flow on arc ( SN[ a ] , EN[ a ] ).
  *
  * The latter OneVarConstraint have fixed 0 LHS and a generic RHS, which can
  * be Inf<Fnumber>(). If *all* the RHS are +Infty, it is possible to use a
  * std::vector<NNConstraint> to represent them instead of a
  * std::vector<LB0Constraint>. The parameter stcc is used to decide if this
  * is done: if
  *
  * - all the RHS are +Infty;
  *
  * - either stcc is not nullptr and it is a SimpleConfiguration<int>;
  *
  * - or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_static_constraints_Configuration is not nullptr,
  *   and it is a SimpleConfiguration<int>;
  *
  * - the f_value of the SimpleConfiguration<int> is != 0
  *
  * then NNConstraint are used to implement *all* bound constraint. Note that
  * this *makes it impossible to change any RHS*, as this would require
  * changing the static Constraint and this is not allowed. Indeed,
  * NNConstraint throws exception if one tries to change its RHS (and LHS as
  * well, but this also NNConstraint does). Hence, if the abstract Constraint
  * are constructed, changing the RHS is not allowed.
  *
  * Note that changing *all* the other parts of any of the FRowConstraint,
  * such as the coefficients of the LinearFunction inside, is not allowed:
  * the MCFBlock will throw exception while processing the corresponding
  * "abstract" Modification. */
 
 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
  override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the MCF
 /** Method that generates the objective of the MCF. Although this would seem
  * to be an exceedingly simple object, there is still a nontrivial decision
  * to be made about it, i.e., whether it is represented as a "sparse"
  * LinearFunction or a "dense" one. This is governed by objc: if
  *
  * - either objc is not nullptr and it is a SimpleConfiguration<double>;
  *
  * - or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_objective_Configuration is not nullptr,
  *   and it is a SimpleConfiguration<double>;
  *
  * then the f_value of the SimpleConfiguration<int> is taken as the "sparsity
  * parameter" (sprs) of the objective function; otherwise, sprs == 0. The
  * parameter is used as follows: if the number of nonzero arc cost
  * coefficients (at the time in which the method is called) is >=
  * sprs * get_NArcs(), then the LinearFunction in the objective is created
  * "dense": each flow Variable is "active" in it, even if it has a zero cost
  * coefficient. Otherwise, the LinearFunction in the objective is created
  * "sparse": only flow Variable with nonzero coefficient are "active" in it.
  * A "dense" Objective makes it much easier to change the cost coefficients
  * (see chg_cost[s]()), but it comes at the cost of more memory. Besides,
  * Solver using it and "trusting" the MCFBlock about how many nonzeroes are
  * there in the LinearFunction may be sorely disappointed, which may have
  * adverse effects on efficiency.
  *
  * Yet, a "dense" Objective is the default, as with sprs == 0 the Objective
  * is created "dense" even if all arc cost coefficients are zero.
  *
  * Note that the decision is taken at the moment in which this method is
  * called, and never changed later, even if the number of nonzeroes
  * changes dramatically. Also, note that if all cost coefficients are
  * "naturally" nonzero, then the Objective will be "dense" no matter what the
  * value of sprs is. Although this may seem obvious, this also means that the
  * Objective will remain "dense" even if later on many coefficients become
  * zero.
  *
  * IMPORTANT NOTE: ALLOWING SPARSE Objective MAKES IT INORDINATELY MORE
  * DIFFICULT TO REACT TO ABSTRACT Modification, WHILE ITS ACTUAL IMPACT ON
  * PERFORMANCES IS VERY DOUBIOUS. THEREFORE, THE SUPPORT FOR IT IS ONLY
  * HALF-BAKED, AND WHATEVER THERE IS IS CURRENTLY COMMENTED OUT. DEVELOPMENT
  * OF THIS FEATURE WILL ONLY BE RESUMED IF CLEAR PROOF OF ITS WORTHINESS
  * IS ACHIEVED.
  *
  * The consequence is that, currently, THE ONLY Modification POSSIBLE TO THE
  * Objective ARE CHANGING THE COEFFICIENTS: DELETING Variable (AND,
  * THEREFORE, ADDING THEM) IS NOT ALLOWED, the MCFBlock will throw exception
  * while processing the corresponding "abstract" Modification. */

 virtual void generate_objective( Configuration *objc = nullptr )
  override final;

/*@} -----------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the MCFBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the MCFBlock
 *  @{ */

 /// get the number of nodes
 inline Index get_NNodes( void ) const { return( NNodes ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of arcs
 inline Index get_NArcs( void ) const { return( SN.size() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of starting nodes
 inline c_Vec_Index & get_SN( void ) const { return( SN ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the starting node of arc i (0 <= i < get_NArcs())
 inline Index get_SN( c_Index i ) const { return( SN[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of ending nodes
 inline  c_Vec_Index & get_EN( void ) const { return( EN ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the ending node of arc i (0 <= i < get_NArcs())
 inline Index get_EN( c_Index i ) const { return( EN[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of arc costs
 /** Returns a const reference to the vector of arc costs. Note that the
  * returned vector can either be of size get_NArcs() or be of size 0, in
  * which case all arc costs are assumed to be 0. */

 inline c_Vec_CNumber & get_C( void ) const { return( C ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the cost of arc i (0 <= i < get_NArcs())
 inline CNumber get_C( c_Index i ) const { return( C.size() ? C[ i ] : 0 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of arc upper bounds
 /** Returns a const reference to the vector of arc upper bounds. Note that
  * the returned vector can either be of size get_NArcs() or be of size 0, in
  * which case all arc upper bounds are assumed to be +Inf. */

 inline c_Vec_FNumber & get_U( void ) const { return( U ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the upper bound of arc i (0 <= i < get_NArcs())

 inline FNumber get_U( c_Index i ) const { return( U.size() ? U[ i ] :
						   Inf<FNumber>() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of node deficits
 /** Returns a const reference to the vector of node deficits. Note that the
  * returned vector can either be of size get_NNodes() or be of size 0, in
  * which case all node deficits are assumed to be 0. Also, note that the
  * position i (0 <= i < get_NNodes()) in this vector correspond to the node
  * whose name is i + 1 as returned from get_SN() and get_EN(). */

 inline c_Vec_FNumber & get_B( void ) const { return( B ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the upper deficit of node i (0 <= i < get_NNodes())
 /** Returns the deficit of node i. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes(). */

 inline FNumber get_B( c_Index i ) const { return( B.size() ? B[ i ] : 0 ); }

/*@} -----------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
    @{ */

 /// returns true if the current solution is (approximately) flow feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the MCFBlock is approximately feasible w.r.t. the flow
  * conservation constraints only. This clearly requires the Variable of the
  * MCFBlock to have been defined, i.e., that generate_abstract_variables() has
  * been called prior to this method. The parameter feps is the relative
  * accuracy defining "approximately". The parameter "useabstract" has the
  * same meaning as in is_feasible() and is_optimal(). Of course, if
  * useabstract == true, then the Constraint of the MCFBlock need to have
  * been generated by calling generate_abstract_constraints() prior to this
  * method. */

 bool flow_feasible( c_FNumber feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) bound feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the MCFBlock is approximately feasible w.r.t. the bound
  * constraints only. This clearly requires the Variable of the MCFBlock to
  * have been defined, i.e., that generate_abstract_variables() has been
  * called prior to this method. The parameter feps is the relative accuracy
  * defining "approximately". The parameter "useabstract" has the same
  * meaning as in is_feasible() and is_optimal(). Of course, if useabstract
  * == true, then the Constraint of the MCFBlock need to have been generated
  * by calling generate_abstract_constraints() prior to this method. */

 bool bound_feasible( c_FNumber feps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) dual feasible
 /** Returns true if the dual solution encoded in the current value of the
  * dual multipliers of both the flow conservation and bound constraints is
  * feasible. This clearly requires the Constraint of the MCFBlock to have
  * been generated by calling generate_abstract_constraints() prior to this
  * method. The parameter ceps is the relative accuracy defining
  * "approximately". The parameter "useabstract" has the same meaning as in
  * is_feasible() and is_optimal(). */

 bool dual_feasible( c_CNumber ceps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if complementary slackness are (approximately) satisfied
 /** Returns true if the (primal) solution encoded in the current value of
  * the flow (x) Variable of the MCFBlock and the dual solution encoded in
  * the current value of the dual multipliers of both the bound constraints
  * approximatively satisfy the Complementary Slackness Conditions. This
  * clearly requires both the Variable and the Constraint of the MCFBlock to
  * have been generated by calling generate_abstract_variables() and
  * generate_abstract_constraints() prior to this method. The parameters ceps
  * and feps are the relative accuracy defining "approximately" respectively
  * for "the reduced cost is zero" and "the flow is at the upper/lower bound".
  * The parameter "useabstract" has the same meaning as in is_feasible() and
  * is_optimal(). */

 bool complementary_slackness( c_CNumber ceps , c_FNumber feps ,
			       bool useabstract = false );

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the MCFBlock is approximately feasible. This clearly
  * requires the Variable of the MCFBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method. Also,
  * if useabstract == true, then the Constraint of the MCFBlock need to have
  * been generated by calling generate_abstract_constraints() prior to this
  * method.
  *
  * The parameter for deciding what "approximately feasible" exactly means is
  * a single FNumber value, representing the *relative* tolerance for
  * satisfaction of both flow conservation constraint and flow upper/lower
  * bounds. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration<FNumber>, then it
  *   if fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration<FNumber>, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0. */
 
 virtual bool is_feasible( bool useabstract = false ,
			   Configuration *fsbc = nullptr ) override final;

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) optimal
 /** Returns true if the solution encoded in the current value of the flow
  * (x) Variable of the MCFBlock is approximately optimal, which means that
  * it is approximately feasible, that the dual solution encoded in the
  * current value of the dual multipliers of both the flow conservation and
  * bound constraints is approximately feasible, and that the two
  * approximately satisfies the Complementary Slackness Conditions. This
  * clearly requires that both the Variable and the Constraint of the
  * MCFBlock to have been defined, i.e., that generate_abstract_variables() and
  * generate_abstract_constraints() have been called prior to this method.
  *
  * This requires two parameters for deciding what "approximately feasible"
  * means, one for the primal (feps) and one for the dual (ceps), like in
  * complementary_slackness(). These are found as follows:
  *
  * - if optc is not nullptr and it is a 
  *   SimpleConfiguration< std::pair<CNumber,FNumber> >, then
  *   ceps = optc->f_value.first and feps = optc->f_value.second;
  *
  * - if optc is not nullptr and it is a SimpleConfiguration<CNumber>, then
  *   ceps = optc->f_value, while feps is taken out of
  *   f_BlockConfig->f_is_feasible_Configuration as in is_feasible();
  *
  * - otherwise, if f_BlockConfig is not nullptr, then feps is taken
  *   out of f_BlockConfig->f_is_feasible_Configuration, while ceps
  *   is taken out of f_BlockConfig->f_is_optimal_Configuration
  *   assuming the latter is a SimpleConfiguration<CNumber>;
  *
  * - otherwise, ceps == feps == 0. */
 
 virtual bool is_optimal( bool useabstract = false  ,
			  Configuration *optc = nullptr ) override final;

/*@} -----------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for R3 Blocks
    @{ */

 /// gets an R3 Block of MCFBlock currently only the copy one
 /** Gets an R3 Block of the MCFBlock. The list of currently supported R3
  * Block is:
  *
  * - r3bc == nullptr: the copy (an MCFBlock identical to the current one).
  *
  */

 virtual Block * get_R3_Block( Configuration *r3bc = nullptr ) override final;

/*--------------------------------------------------------------------------*/
 /// maps back the solution from a copy MCFBlock to the current one
 /** Maps back the solution from a copy MCFBlock to the current one. The
  * parameter r3bc is useless (has to be nullptr). The parameter solc decides
  * which part of the solution is mapped:
  *
  * - if solc != nullptr and it is a SimpleConfiguration<int>, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration != nullptr and it
  *   is a SimpleConfiguration<int>, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that R3B may not contain some or all of the required solution, if
  * the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception. */ 

 virtual void map_back_solution( Block *R3B , Configuration *r3bc = nullptr ,
				 Configuration *solc = nullptr )
  override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps the solution of the current MCFBlock to a copy MCFBlock
 /** Maps the solution of the current MCFBlock to a copy MCFBlock. The
  * parameter r3bc is useless (has to be nullptr). The parameter solc decides
  * which part of the solution is mapped:
  *
  * - if solc != nullptr and it is a SimpleConfiguration<int>, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration != nullptr and it
  *   is a SimpleConfiguration<int>, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that the current MCFBlock may not contain some or all of the
  * required solution, if the corresponding Variable/Constraint have not
  * been constructed yet: this throws an exception. */ 

 virtual void map_forward_solution( Block *R3B ,
				    Configuration *r3bc = nullptr ,
				    Configuration *solc = nullptr )
  override final;

/*--------------------------------------------------------------------------*/
 /** No specific Configuration is required, hence expected, for MCFBlock.
  *
  * IMPORTANT NOTE: map_forward_Modification() only maps "physical"
  * Modification. The point is that if any part of the "abstract
  * representation" of MCFBlock is changed, the corresponding "abstract"
  * Modification is intercepted in add_Modification() and a "physical"
  * Modification is also issued. Hence, for any change in MCFBlock there
  * will always be both Modification "in flight", and therefore there is
  * no need (and good reasons not) to map both.
  *
  * In particular, the method handles the following Modification:
  *
  * - GroupModification
  *
  * - MCFBlockRngdMod
  *
  * - MCFBlockSbstMod
  *
  * - NBModification
  *
  * Any other Modification is ignored (and false is returned).
  *
  * Note that for GroupModification, true is returned only if all the
  * inner Modification of the GroupModification return true.
  *
  * Note that if the issueAMod param is eModBlck, then it is "downgraded" to
  * eNoBlck: the method directly does "physical" changes, hence there is no
  * reason for it to issue "abstract" Modification with concerns_Block() ==
  * true. */

 virtual bool map_forward_Modification( Block *R3B , sp_Mod mod ,
					Configuration *r3bc = nullptr ,
					c_ModParam issuePMod = eNoBlck ,
					c_ModParam issueAMod = eModBlck )
  override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** No specific Configuration is required, hence expected, for MCFBlock.
  *
  * The current implementation of map_back_Modification() actually uses
  * map_forward_Modification() in reverse, so see the comments to the latter
  * method. */

 virtual bool map_back_Modification( Block *R3B , sp_Mod mod ,
				     Configuration *r3bc = nullptr ,
				     c_ModParam issuePMod = eNoBlck ,
				     c_ModParam issueAMod = eModBlck )
  override final;

/*@} -----------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a MCFSolution representing the current solution of this MCFBlock
 /** Returns a MCFSolution representing the current solution status of this
  * MCFBlock.The parameter solc decides which part of the solution is saved:
  *
  * - if solc != nullptr and it is a SimpleConfiguration<int>, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration != nullptr and it
  *   is a SimpleConfiguration<int>, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that MCFBlock may not contain some or all of the required solution,
  * if the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception, unless emptys = true, in which case thew
  * MCFSolution object is only prepped for getting a solution, but it is not
  * really getting one now.
  *
  * Note that, although the method clearly returns a MCFSolution, formally
  * the return type is Solution *. This is because it is not possible to
  * forward declare MCFSolution as a derived class from Solution, nor to
  * define MCFSolution before MCFBlock because the former uses some type
  * information declared in the latter. */ 

 virtual Solution * get_Solution( Configuration *solc = nullptr ,
				  bool emptys = true ) override final;

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the flow solution
 /** Method to get the flow solution; upon return, FSol[ i ] contains the
  * current value of the flow solution for arc strt + i for all 0 <= i < 
  * min( stp , get_NArcs() ). */

 void get_x( Vec_FNumber & FSol , c_Index strt = 0 ,
	                          c_Index stp = Inf<Index>() );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the flow solution for an arbitrary subset of arcs
 /** Method to get the flow solution; upon return, FSol[ i ] contains the
  * current value of the flow solution for arc nms[ i ] for all 0 <= i < 
  * nms.size(). */

 void get_x( Vec_FNumber & FSol , c_Vec_Index & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the flow solution of the given arc

 FNumber get_x( c_Index arc ) {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  return( x[ arc ].get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the potential solution
 /** Method to get the potential solution; upon return, PSol[ i ] contains the
  * current value of the potential solution for node strt + i for all 0 <= i <
  * min( stp , get_NNodes() ). Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes(). */

 void get_pi( Vec_CNumber & PSol , c_Index strt = 0 ,
	                           c_Index stp = Inf<Index>() );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the flow potential for an arbitrary subset of nodes
 /** Method to get the potential solution; upon return, PSol[ i ] contains the
  * current value of the potential solution for node nms[ i ] for all 0 <= i
  * < nms.size(). Note that "node names" here go from 0 to get_NNodes() - 1,
  * despite the fact that get_SN() and get_EN() report node "names" between
  * 1 and get_NNodes(). */

 void get_pi( Vec_CNumber & PSol , c_Vec_Index & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the potential solution of the given node
 /** Method to get the potential solution of the given node; note that "node
  * names" here go from 0 to get_NNodes() - 1, despite the fact that get_SN()
  * and get_EN() report node "names" between 1 and get_NNodes(). */

 CNumber get_pi( c_Index nde ) {
  if( ! E.size() )
   throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );
  if( nde >= get_NNodes() )
   throw( std::invalid_argument( "invalid node name" ) );
   
  return( E[ nde ].get_dual() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the reduced costs
 /** Method to get the reduced costs; upon return, RC[ i ] contains the
  * current value of the reduced cost for arc strt + i for all 0 <= i < 
  * min( stp , get_NArcs() ). */

 void get_rc( Vec_CNumber & RC , c_Index strt = 0 ,
	                         c_Index stp = Inf<Index>() );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the reduced costs for an arbitrary subset of arcs
 /** Method to get the reduced costs; upon return, RC[ i ] contains the
  * current value of the reduced costs for arc nms[ i ] for all 0 <= i < 
  * nms.size(). */

 void get_rc( Vec_CNumber & RC , c_Vec_Index & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the reduced costs of the given arc

 CNumber get_rc( c_Index arc ) {
  if( ! E.size() )
   throw( std::logic_error( "reduced costs unavailable if Constraint aren't"
			    ) );
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  auto nnc = boost::any_cast<std::vector<NNConstraint> *>(
					     get_static_constraints()[ 1 ] );
  if( nnc )
   return( (*nnc)[ arc ].get_dual() );

  auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					     get_static_constraints()[ 1 ] );
  assert( lbc );
  return( (*lbc)[ arc ].get_dual() );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the flow solution
 /** Method to set the flow solution; the values found in the c_Vec_FNumber
  * between fstrt (included) and fstop (excluded) are copied into the value of
  * the flow variable x[ strt + i ]. This is typically used by a Solver. */

 void set_x( c_Vec_FNumber_it fstrt , c_Vec_FNumber_it fstop ,
	     c_Index strt = 0 )
 {
  if( std::distance( fstop , fstrt ) + strt > get_NArcs() )
   throw( std::invalid_argument( "too many values provided" ) );

  for( auto xi = x.begin() + strt ; fstrt < fstop ; )
   (xi++)->set_value( *(fstrt++) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the flow solution of the given arc

 void set_x( c_Index arc , c_FNumber FSol ) {
  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  x[ arc ].set_value( FSol );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the potential solution
 /** Method to set the potential solution; the values found in the
  * c_Vec_CNumber between pstrt (included) and pstop (excluded) are copied
  * into the potential of node (dual multiplier of the flow balance
  * constraint) strt + i. This is typically used by a Solver. */

 void set_pi( c_Vec_CNumber_it pstrt , c_Vec_CNumber_it pstop ,
	      c_Index strt = 0 )
 {
  if( ! E.size() )  // nowhere to put the value
   return;          // cowardly (and silently) return

  if( std::distance( pstop , pstrt ) + strt > get_NNodes() )
   throw( std::invalid_argument( "too many values provided" ) );

  for( auto Ei = E.begin() + strt ; pstrt < pstop ; )
   (Ei++)->set_dual( *(pstrt++) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the potential solution of the given node

 void set_pi( CNumber PSol , c_Index nde ) {
  if( ! E.size() )  // nowhere to put the value
   return;          // cowardly (and silently) return

  if( nde >= get_NNodes() )
   throw( std::invalid_argument( "invalid node name" ) );

  return( E[ nde ].set_dual( PSol ) );
  }

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the reduced costs
 /** Method to set the reduced costs solution; the values found in the
  * c_Vec_CNumber between rcstrt (included) and rcstop (excluded) are copied
  * into the reduced cost of arc (dual value of the bound constraint) strt +
  * i. This is typically used by a Solver. */

 void set_rc( c_Vec_CNumber_it rcstrt , c_Vec_CNumber_it rcstop ,
	      c_Index strt = 0 )
 {
  if( ! E.size() )  // nowhere to put the value
   return;          // cowardly (and silently) return

  if( std::distance( rcstop , rcstrt ) + strt > get_NArcs() )
   throw( std::invalid_argument( "too many values provided" ) );

  auto nnc = boost::any_cast<std::vector<NNConstraint> *>(
					     get_static_constraints()[ 1 ] );
  if( nnc ) {
   for( auto bi = nnc->begin() + strt ; rcstrt < rcstop ; )
    (bi++)->set_dual( *(rcstrt++) );
   }
  else {
   auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					     get_static_constraints()[ 1 ] );
   assert( lbc );
   for( auto bi = lbc->begin() + strt ; rcstrt < rcstop ; )
    (bi++)->set_dual( *(rcstrt++) );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the reduced cost of the given arc

 void set_rc( c_CNumber RC , c_Index arc ) {
  if( ! E.size() )  // nowhere to put the value
   return;          // cowardly (and silently) return

  if( arc >= get_NArcs() )
   throw( std::invalid_argument( "invalid arc name" ) );

  auto nnc = boost::any_cast<std::vector<NNConstraint> *>(
					     get_static_constraints()[ 1 ] );
  if( nnc )
   (*nnc)[ arc ].set_dual( RC );
  else {
   auto lbc = boost::any_cast<std::vector<LB0Constraint> *>(
					     get_static_constraints()[ 1 ] );
   assert( lbc );
   (*lbc)[ arc ].set_dual( RC );
   }
  }

/*@} -----------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the MCFBlock
 /** Method for handling Modification.
  *
  * The version of MCFBlock has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the MCFBlock, and "translate"
  * them into both changes of the actual data structures and corresponding
  * "physical Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the MCFBlock, and this must not be
  * repeated. In particular, this would be an issue if the Modification would
  * be [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * MCFBlock), which would mean that the "physical Modification" would be
  * issued twice.
  *
  * The following "abstract Modification" are handled:
  *
  * - GroupModification, that are simply unpacked into the individual
  *   sub-[Group]Modification and dealt with individually;
  *
  * - LinearFunctionModSbst adding/removing Variable and changing coefficients
  *   coming from the (LinearFunction into the FRow)Objective, but *not* from
  *   the (LinearFunction into the FRow)Constraint;
  *
  * - LinearFunctionModRngd adding/removing Variable and changing coefficients
  *   coming from the (LinearFunction into the FRow)Objective, but *not* from
  *   the (LinearFunction into the FRow)Constraint;
  *
  * - RowConstraintMod changing the RHS of the bound constraints and both
  *   sides at once of the flow conservation ones, but not any other
  *   combination (and note that the RHS of the bound constraints may not
  *   be changeable at all depending on how they have been constructed; in
  *   this case attempting to do this throws exception, which means that there
  *   is no need to handling this here);
  *
  * - VariableMod fixing and un-fixing a flow ColVariable; however, note
  *   that *fixing is only permitted if the value() of the ColVvariable is
  *   zero*, because that corresponds to closing the arc, exception being
  *   thrown otherwise.
  *
  * Any other Modification reaching the MCFBlock will lead to exception
  * being thrown. */

 virtual void add_Modification( sp_Mod mod , ChnlName chnl = 0 )
  override final;

/*@} -----------------------------------------------------------------------*/
/*---------------------- Methods for handling Solver -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solver
    @{ */


/*@} -----------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE MCFBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the MCFBlock
 *  @{ */

 virtual void serialize( netCDF::NcGroup&& group ) const override final;

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the MCF instabnce
 *
 * All the methods in this section have two parameters issueMod and issueAMod
 * which control if and how the, respectively, "physical Modification" and
 * "abstract Modification" corresponding to the change have to be issued, and
 * where (to which channel). The format of the parameters is that of
 * Observer::make_par(), except that the value eModBlck is ignored and
 * treated it as if it were eNoBlck [see Observer::issue_pmod()]. This is
 * because it makes no sense to issue an "abstract" Modification with
 * concerns_Block() == true, since the changes in the MCFBlock have surely
 * been done already, and this is just not possible for a "physical"
 * Modification.
 *
 * IMPORTANT NOTE: the current implementation of all these methods issues (at
 * most) *two separate* Modification, a "physical" and an "abstract" one. The
 * latter may be a GroupModification bunching together related abstract
 * Modification, but the two Modification are nonetheless separate. A
 * different approach could be to issue a single GroupModification with inside
 * both the "physical" and the "abstract" one (the latter possibly itself a
 * GroupModification). This may allow a more efficient handling of
 * Modification by ensuring that the two are always received together, but at
 * the cost of a more intricate code that is best avoided for now.
 *
 * Note: the methods accept the eDryRun value for the issueAMod parameter for
 * the "abstract" representation. This allows to re-use them within MCFBlock
 * itself when reacting to abstract Modification, where the  "abstract"
 * representation has been changed already. However, the eDryRun value is not
 * allowed (it is ignored) for the issuePMod parameter for the "physical"
 * representation, as there is no reasonable use for this. Basically, this
 * makes eDryRun equivalent to eNoMod.
 *  @{ */

 /// change the costs of a contiguous interval of arcs
 /** Method to change the costs of a subset of arcs with "contiguous names".
  * That is, *( NCost + i - strt ) becomes the new cost of arc i for all
  * strt <= i < min( stop , get_NArcs() ).
  *
  * Note that, if the Objective is a "sparse" LinearFunction (see
  * compute_objective()), then changing the costs can issue up to three
  * different Modification; in particular a LinearFunctionMod for adding a
  * Variable (setting to nonzero a previously zero coefficient), one for
  * removing Variable (vice-versa), and one for modifying the coefficients.
  * If more than one Modification is actually issued and issueAMod specifies
  * an open channel, then the channel is nested so that the three Modification
  * are grouped into a single GroupModification. Similarly, if instead
  * issueAMod specifies the default channel, then a new channel is opened to
  * group the multiple Modification and immediately closed when the last one
  * is issued. If, instead, the Objective is a "dense" LinearFunction, then
  * at most one LinearFunctionMod for modifying the coefficients is issued.
  * Of course this only applies if issueAMod specifies that abstract
  * Modification have to be issued *and* the abstract Objective has been
  * constructed.
  *
  * Also, if issueMod says so then a "physical" MCFBlockRngdMod is issued. */

 void chg_costs( c_Vec_CNumber_it NCost , c_Index strt = 0 ,
		 Index stop = Inf<Index>() , c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the costs of an arbitrary subset of arcs
 /** Method to change the costs of an arbitrary subset of arc. That is,
  * *( NCost + i ) becomes the new cost of arc nms[ i ] for all 0 <= i <
  * NCost.size(), (which means that nms.size() == NCost.size()). The
  * parameter ordered tells if the nms vector is ordered for increasing
  * index of the arc. As the the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate MCFBlockSbstMod object.
  *
  * See chg_costs( range ) for Modification issued (except that, of course,
  * the "physical" one is a MCFBlockSbstMod). */

 void chg_costs( c_Vec_CNumber_it NCost , Vec_Index && nms ,
		 const bool ordered = false ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the cost of the given arc
 /** Changes the cost of the given arc.
  *
  * Note that this can issue only one Modification of each type; the
  * "physical" one is a MCFBlockRngdMod with stop = start + 1. */

 void chg_cost( c_CNumber NCost , c_Index arc , 
		c_ModParam issueMod = eNoBlck ,
		c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the capacities of a contiguous interval of arcs
 /** Method to change the capacities of a subset of arcs with "contiguous
  * names". That is, *( NCap + i - strt ) becomes the new capacity of arc i
  * for all strt <= i < min( stop , get_NArcs() ). Note that, according to
  * the Configuration of the static Constraint, the capacity of the arcs
  * cannot be changed: trying to do that will result in an exception being
  * thrown.
  *
  * Note that changing the capacities can issue as many Modification as there
  * are arcs in the range, in particular OneVarConstraintMod with type
  * RowConstraintMod::eChgRHS. If more than one Modification is actually
  * issued and issueAMod specifies an open channel, then the channel is
  * nested so that all the Modification are grouped into a single
  * GroupModification. Similarly, if instead issueAMod specifies the default
  * channel, then a new channel is opened to group the multiple Modification
  * and immediately closed when the last one is issued. Of course this only
  * applies if issueAMod specifies that abstract Modification have to be
  * issued *and* the abstract Constraint have been constructed.
  *
  * Also, if issueMod says so then a "physical" MCFBlockRngdMod is issued. */

 void chg_ucaps( c_Vec_FNumber_it NCap , c_Index strt = 0 ,
		 Index stop = Inf<Index>() , c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the capacities of an arbitrary subset of arcs
 /** Method to change the capacities of an arbitrary subset of arc. That is,
  * *( NCap + i ) becomes the new capacity of arc nms[ i ] for all 0 <= i <
  * NCap.size() (which means that nms.size() == NCap.size()). The parameter
  * ordered tells if the nms vector is ordered for increasing index of the
  * arc. As the the && tells, nms is "consumed" by the method, typically
  * being shipped to an appropriate MCFBlockSbstMod object.
  *
  * Note that, according to the Configuration of the static Constraint, the
  * capacity of the arcs cannot be changed: trying to do that will result in
  * an exception being thrown.
  *
  * See chg_ucaps( range ) for Modification issued (except that, of course,
  * the "physical" one is a MCFBlockSbstMod). */

 void chg_ucaps( c_Vec_FNumber_it NCap , Vec_Index && nms ,
		 const bool ordered = false ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the capacity of the given arc
 /** Method to change the capacity of a given arc: NCap becomes the new
  * capacity of arc arc. Note that, according to the Configuration of the
  * static Constraint, the capacity of the arcs cannot be changed: trying to
  * do that will result in an exception being thrown.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * MCFBlockRngdMod with stop = start + 1. */

 void chg_ucap( c_FNumber NCap , c_Index arc ,
		c_ModParam issueMod = eNoBlck ,
		c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the deficits of a contiguous interval of nodes
 /** Method to change the deficits of a subset of nodes with "contiguous
  * names". That is, *( NDfct + i - strt ) becomes the new deficit of node i
  * for all strt <= i < min( stop , get_NNodes() ). Note that "node names"
  * here go from 0 to get_NNodes() - 1, despite the fact that get_SN() and
  * get_EN() report node "names" between 1 and get_NNodes().
  *
  * Note that changing the capacities can issue as many Modification as there
  * are nodes in the range, in particular FRowConstraintMod with type
  * RowConstraintMod::eChgBTS. If more than one Modification is actually
  * issued and issueAMod specifies an open channel, then the channel is
  * nested so that all the Modification are grouped into a single
  * GroupModification. Similarly, if instead issueAMod specifies the default
  * channel, then a new channel is opened to group the multiple Modification
  * and immediately closed when the last one is issued. Of course this only
  * applies if issueAMod specifies that abstract Modification have to be
  * issued *and* the abstract Constraint have been constructed.
  *
  * Also, if issueMod says so then a "physical" MCFBlockRngdMod is issued. */

 void chg_dfcts( c_Vec_FNumber_it NDfct , c_Index strt = 0 ,
		 Index stop = Inf<Index>() ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the deficits of an arbitrary subset of nodes
 /** Method to change the deficits of an arbitrary subset of nodes. That is,
  * *( NDfct + i ) becomes the new deficit of node nms[ i ] for all 0 <= i <
  * NDfct.size(), (which means that nms.size() == NDfct.size()). The
  * parameter ordered tells if the nms vector is ordered for increasing index
  * of the node. Note that "node names" here go from 0 to get_NNodes() - 1,
  * despite the fact that get_SN() and get_EN() report node "names" between
  * 1 and get_NNodes(). As the the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate MCFBlockSbstMod object.
  *
  * See chg_dfcts( range ) for Modification issued (except that, of course,
  * the "physical" one is a MCFBlockSbstMod). */

 void chg_dfcts( c_Vec_FNumber_it NDfct , Vec_Index && nms ,
		 const bool ordered = false ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the deficit of the given node
 /** Method to change the deficit of a given node: NDfct becomes the new
  * deficit of node nde. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_SN() and get_EN() report node
  * "names" between 1 and get_NNodes().
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * MCFBlockRngdMod with stop = start + 1. */

 void chg_dfct( c_FNumber NDfct , c_Index nde ,
		c_ModParam issueMod = eNoBlck ,
		c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// closes a contiguous interval of arcs
 /** Method to close a subset of arcs with all "contiguous names" comprised
  * between strt (included) and min( stop , get_NArcs() ) (excluded). The 
  * flow on the arcs is fixed to 0 but the arcs are not removed from the
  * problem, and their capacity and cost are not changed, so that they can be
  * easily re-opened later. When the problem is created, all arcs are open.
  * Closing an already closed arc does nothing.
  *
  * Note that closing multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod with type
  * Variable::kFixed. If more than one Modification is actually issued and
  * issueAMod specifies an open channel, then the channel is nested so that
  * all the Modification are grouped into a single GroupModification.
  * Similarly, if instead issueAMod specifies the default channel, then a
  * new channel is opened to group the multiple Modification and immediately
  * closed when the last one is issued. Of course this only applies if
  * issueAMod specifies that abstract Modification have to be issued.
  *
  * Also, if issueMod says so then a "physical" MCFBlockRngdMod is issued. */

 void close_arcs( c_Index strt = 0 , Index stop = Inf<Index>() ,
		  c_ModParam issueMod = eNoBlck ,
		  c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// closes an arbitrary subset of arcs
 /** Method to close an arbitrary subset of arc, i.e., all those whose names
  * are found in the array nms. The flow on the arcs is fixed to 0 but the
  * arcs are not removed from the problem, and their capacity and cost are
  * not changed, so that they can be easily re-opened later. When the problem
  * is created, all arcs are open. Closing an already closed arc does
  * nothing.
  *
  * The parameter ordered tells if the nms vector is ordered for increasing 
  * index of the arc. As the the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate MCFBlockSbstMod object.
  *
  * See close_arcs( range ) for Modification issued (except that, of course,
  * the "physical" one is a MCFBlockSbstMod). */

 void close_arcs( Vec_Index && nms , const bool ordered = false ,
		  c_ModParam issueMod = eNoBlck ,
		  c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// closes the given arc
 /** Method to "close" the given arc: the flow on arc is fixed to 0. The arc
  * is not removed from the problem, and its capacity and cost are not
  * changed, so that it can be easily re-opened later. When the problem is
  * created, all arcs are open. Closing an already closed arc does nothing.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * MCFBlockRngdMod with stop = start + 1. */

 void close_arc( c_Index arc ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/// re-opens a contiguous interval of arcs
 /** Method to "open" a subset of closed arcs with all "contiguous names"
  * comprised between strt (included) and min( stop , get_NArcs() )
  * (excluded). Opening an already open arc (which is what all arcs are when
  * the problem is created) does nothing.
  *
  * Note that opening multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod with type
  * ColVariable::kContinuous. If more than one Modification is actually
  * issued and issueAMod specifies an open channel, then the channel is
  * nested so that all the Modification are grouped into a single
  * GroupModification. Similarly, if instead issueAMod specifies the default
  * channel, then a new channel is opened to group the multiple Modification
  * and immediately closed when the last one is issued. Of course this only
  * applies if issueAMod specifies that abstract Modification.
  *
  * Also, if issueMod says so then a "physical" MCFBlockRngdMod is issued. */

 void open_arcs( c_Index strt = 0 , Index stop = Inf<Index>() ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// re-opens an arbitrary subset of arcs
 /** Method to "open" an arbitrary subset of closed arc, i.e., all those
  * whose names are found in the array nms. Opening an already open arc
  * (which is what all arcs are when the problem is created) does nothing.
  *
  * The parameter ordered tells if the nms vector is ordered for increasing 
  * index of the arc. As the the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate MCFBlockSbstMod object.
  *
  * Note that closing multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod with type
  * Variable::kFixed. If more than one Modification is actually issued and
  * issueAMod specifies an open channel, then the channel is nested so that
  * all the Modification are grouped into a single GroupModification.
  * Similarly, if instead issueAMod specifies the default channel, then a
  * new channel is opened to group the multiple Modification and immediately
  * closed when the last one is issued.
  * Of course this only applies if issueAMod specifies that abstract
  * Modification have to be issued *and* the abstract Constraint have been
  * constructed.
  *
  * Also, if issueMod says so then a "physical" MCFBlockRngdMod is issued. */

 void open_arcs( Vec_Index && nms , const bool ordered = false ,
		 c_ModParam issueMod = eNoBlck ,
		 c_ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// re-opens the given arc
 /** Method to "open" the given closed arc, i.e., allow the flow on arc to
  * vary. Opening an already open arc (which is what all arcs are when the
  * problem is created) does nothing.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * MCFBlockRngdMod with stop = start + 1. */

 void open_arc( c_Index arc ,
		c_ModParam issueMod = eNoBlck ,
		c_ModParam issueAMod = eNoBlck );

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED FRIENDS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
    @{ */

 /// print the MCFBlock on an ostream with the given verbosity
 /** Protected method to print information about the MCFBlock; with the
  * "complete" level it outputs the MCFBlock in DIMACS formar. */

 virtual void print( std::ostream &output ) const override final;

/*--------------------------------------------------------------------------*/
 /// loads the MCF instance from file in DIMACS standard format
 /** Protected method for loading a MCFBlock out of a std::istream (which is
  * what operator>> is dispatched to. The std::istream is assumed to contain
  * the description of a MCF instance in DIMACS standard format, which is 
  * the following. The first line must be
  *
  *      p min <number of nodes> <number of arcs>
  *
  * Then the node definition lines must be found, in the form
  *
  *      n <node number> <node supply>
  *
  * Not all nodes need have a node definition line; these are given zero
  * supply, i.e., they are transhipment nodes (supplies are the inverse of
  * deficits, i.e., a node with positive supply is a source node). Finally,
  * the arc definition lines must be found, in the form
  *
  *    a <start node> <end node> <lower bound> <upper bound> <flow cost>
  *
  * There must be exactly <number of arcs> arc definition lines in the file.
  *
  * Note that the file format accepted by LoadMCF is more general than the
  * DIMACS standard format, in that node and arc definitions can be mixed in
  * any order, while the DIMACS file requires all node information to appear
  * before all arc information. Also, capacities of arcs can be set to
  * +Inf<FNumber>() by putting "INF", "Inf" or "inf" in the file (actually,
  * any string starting with "I" or "i" where these would be expected.
  *
  * Like load( memory ), if there is any Solver attached to this MCFBlock
  * then a NBModification (the "nuclear option") is issued. */

 virtual void load( std::istream &input ) override final;

/*@}------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 Index NNodes;                        ///< the number of nodes
 Vec_Index SN;                        ///< vector of arc starting nodes
 Vec_Index EN;                        ///< vector of arc ending nodes

 Vec_CNumber C;                       ///< vector of arc costs
 Vec_FNumber U;                       ///< vector of arc upper capacities
 Vec_FNumber B;                       ///< vector of node deficits

 std::vector<ColVariable> x;          ///< the flow variables
 std::vector<FRowConstraint> E;       ///< the flow conservation constraints

 FRealObjective c;                    ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 inline Index p2i( Variable * const var );

 void guts_of_destructor( void );

 void guts_of_add_Modification( sp_Mod mod );

 inline ModParam make_amod_param( c_ModParam issueAMod , c_Index num );

 inline void unmake_amod_param( c_ModParam oldiAM , c_ModParam newiAM ,
				c_Index num );

 inline Index bound_number( Constraint * const Cnst );

 inline Index const_number( FRowConstraint * const Cnst );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;        // insert MCFBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( class( MCFBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFBlockMod -----------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for "simple" modifications to a MCFBlock
/** Derived class from Modification to describe "simple" modifications to a
 *  MCFBlock, which is "everything changed". Note that it is derived from
 *  Modification rather than, say, BlockMod (which has the same structure)
 *  because this is a class of "physical Modification". This means that any
 *  MCFBlockMod refers to changes in the "physical representation" of the
 *  MCFBlock; the corresponding changes in the "abstract representation" of
 *  the MCFBlock are dealt with by means of "abstract Modification", i.e.,
 *  derived classes from AModification (as is BlockMod, which is why
 *  MCFBlockMod cannot derive from BlockMod). */

class MCFBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of MCFBlockMod
 /** Actually, this public enum is not used in the base BlockMod class.
  * However, both MCFBlockRngdMod and MCFBlockSbstMod derive from MCFBlockMod
  * and require it, so it makes sense to define it only once in the base
  * class. */
 
 enum MCFB_mod_type {
  eChgCost = 0 ,   ///< change the arc costs
  eChgCaps     ,   ///< change the arc capacities
  eChgDfct     ,   ///< change the node deficits
  eOpenArc     ,   ///< open arcs
  eCloseArc        ///< close arcs
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 MCFBlockMod( MCFBlock *fblock ) : f_Block( fblock ) {}
 ///< constructor: takes the MCFBlock

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~MCFBlockMod() { }   ///< destructor, does nothing

/*---------------------- PUBLIC FIELDS OF THE CLASS ------------------------*/

 MCFBlock *f_Block;
             ///< pointer to the MCFBlock to which the MCFBlockMod refers

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the MCFBlockMod
 virtual inline void print( std::ostream &output ) const {
  output << "MCFBlockMod[" << this << "]" << std::endl;
  }

/*--------------------------------------------------------------------------*/

 };  // end( class( MCFBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MCFBlockRngdMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from MCFBlockMod for "ranged" modifications
/** Derived class from MCFBlockMod to describe "ranged"
 * modifications to a MCFBlock, i.e., modifications that apply to an interval
 * of either arcs or nodes. */

class MCFBlockRngdMod : public MCFBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the MCFBlock, the type, and the range
 MCFBlockRngdMod( MCFBlock *fblock , int type ,
		  MCFBlock::Index strt , MCFBlock::Index stop )
  : MCFBlockMod( fblock ) , f_type( type ) , f_strt( strt ) , f_stop( stop )
 {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~MCFBlockRngdMod() { }   ///< destructor, does nothing

/*--------------------- PUBLIC FIELDS OF THE CLASS ------------------------*/

 int f_type;                 ///< type of modification

 MCFBlock::Index f_strt;     ///< begin of the range
 MCFBlock::Index f_stop;     ///< end of the range
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the MCFBlockRngdMod
 virtual inline void print( std::ostream &output ) const {
  output << "MCFBlockRngdMod[" << this << "]: ";
  switch( f_type ) {
   case( eChgCost ): output << "change costs "; break;
   case( eChgCaps ): output << "change capacities "; break;
   case( eChgDfct ): output << "change deficits "; break;
   case( eOpenArc ): output << "open arcs "; break;
   default:          output << "close arcs ";
   }
  output << "[ " << f_strt << ", " << f_stop - 1 << " ]" << std::endl;
  }

/*--------------------------------------------------------------------------*/

 };  // end( class( MCFBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MCFBlockSbstMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from MCFBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  MCFBlock, i.e., modifications that apply to an arbitrary subset of either
 * the arcs or the nodes. */

class MCFBlockSbstMod : public MCFBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:


/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the MCFBlock, the type, and the subset
 /**< Constructor: takes the MCFBlock, the type, and the subset. As the the
  * && tells, nms is "consumed" by the constructor and its resources become
  * property of the MCFBlockSbstMod object. */

 MCFBlockSbstMod( MCFBlock *fblock , int type , MCFBlock::Vec_Index && nms )
  : MCFBlockMod( fblock ) , f_type( type ) , f_nms( std::move( nms ) ) { }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~MCFBlockSbstMod() { }   ///< destructor, does nothing

/*--------------------- PUBLIC FIELDS OF THE CLASS ------------------------*/

 int f_type;                  ///< type of modification

 MCFBlock::Vec_Index f_nms;   ///< the subset
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the MCFBlockSbstMod
 virtual inline void print( std::ostream &output ) const {
  output << "MCFBlockSbstMod[" << this << "] ";
  switch( f_type ) {
   case( eChgCost ): output << "change costs "; break;
   case( eChgCaps ): output << "change capacities "; break;
   case( eChgDfct ): output << "change deficits "; break;
   case( eOpenArc ): output << "open arcs "; break;
   default:          output << "close arcs ";
   }
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------------------------------------------------------------*/

 };  // end( class( MCFBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFSolution -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a MCFBlock
/** The MCFSolution class, derived from Solution, represents a solution of a
 * MCFBlock, i.e.:
 *
 * - an m-vector of FNumber for the arc flow values;
 *
 * - an n-vector of CNumber for the node potentials;
 *
 * where m is the number of arcs and n is the number of nodes in the graph.
 * This means that
 *
 *       THE REDUCED COSTS ARE NOT EXPLICITLY SAVED
 *
 * This is OK for feasible dual solutions, as the dual variables of the bound
 * constraints (a.k.a. Reduced Costs) can be cheapily computed out of the
 * potentials. This may not be appropriate in all cases, as one may want to
 * deal with unfeasible dual solutions; if this will ever be the case, the
 * MCFSolution class will have to be changed accordingly.
 *
 * It is useful to remark that some special cases of MCF would actually have
 * "special" solutions ("less general" ones in the parlance of Solution). In
 * particular:
 *
 * - if all capacities are Inf<FNumber>() and there is only one source or sink
 *   node, then the MCF problem is in fact a Shortest Path (sub-)Tree one, and
 *   its solutions can be represented by means of a predecessor function;
 *
 * - if all (finite) capacities and node deficits are integer, then there
 *   always exist optimal flow solutions of MCF that are integer;
 *
 * - if all arc costs are integer, then there always exist optimal potential
 *   solutions of MCF that are integer.
 *
 * Thus, MCFBlock would have scope for different kinds of Solution objects.
 * The currently implemented one is the "most general" one, so that the
 * Solution::scale() and Solution::sum() operations are always possible;
 * specialized Solution for the specific cases are left for future
 * development. */

class MCFSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend MCFBlock;  ///< make MCFBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING MCFSolution ----------------*/

 MCFSolution() { }  /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~MCFSolution() { }  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A MCFSolution -----------*/

 virtual void read( const Block * const block ) override final;

 virtual void write( Block * const block ) override final;

 virtual MCFSolution * scale( double factor ) const override final;

 virtual void sum( const Solution * solution , double multiplier )
  override final;

 virtual MCFSolution * clone( bool empty = false ) const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 virtual void print( std::ostream &output ) const override final {
  output << "MCFSolution [" << this << "]: " << v_x.size() << " flows and "
	 << v_pi.size() << " potentials" << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 MCFBlock::Vec_FNumber v_x;   ///< the arc flows

 MCFBlock::Vec_CNumber v_pi;  ///< the node potentials

/*--------------------------------------------------------------------------*/

 };  // end( class( MCFSolution ) )

/*@}  end( group( MCFBlock_CLASSES ) ) -------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MCFBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/

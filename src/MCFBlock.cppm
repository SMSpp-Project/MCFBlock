/*--------------------------------------------------------------------------*/
/*--------------------------- File MCFBlock.cppm ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.mcfblock, which exports every name that
 * the headers of MCFBlock declare at namespace scope, together with the
 * modules it depends on; a translation unit can then write
 *
 *     import smspp.mcfblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "MCFBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.mcfblock;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // MCFBlock.h
 using SMSpp_di_unipi_it::MCFBlock;
 using SMSpp_di_unipi_it::MCFBlockMod;
 using SMSpp_di_unipi_it::MCFBlockRngdMod;
 using SMSpp_di_unipi_it::MCFBlockSbstMod;
 using SMSpp_di_unipi_it::MCFSolution;
 using SMSpp_di_unipi_it::Vec_MCFBlock;
 using SMSpp_di_unipi_it::Vec_MCFBlock_it;
 using SMSpp_di_unipi_it::c_Vec_MCFBlock;
 using SMSpp_di_unipi_it::c_Vec_MCFBlock_it;
 using SMSpp_di_unipi_it::p_MCFBlock;
}

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFBlock.cppm -------------------------*/
/*--------------------------------------------------------------------------*/

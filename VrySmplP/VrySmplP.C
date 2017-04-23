/*--------------------------------------------------------------------------*/
/*----------------------- File VrySmplP.C ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Solves the Very Simple Problem in {0, 1} optimization, that is, the one
 * without any constraint. However, also implements an efficient procedure
 * for listing all the solutions of the problem in nonincreasing (for a
 * maximization problem) order of the objective function value.
 *
 * \version 1.01
 *
 * \date 24 - 04 - 2004
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Giovanni Rinaldi \n
 *         Istituto di Analisi di Sistemi e Informatica \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * Copyright &copy 2003 - 2004 by Antonio Frangioni, Giovanni Rinaldi
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- MACROS --------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- INCLUDES -------------------------------------*/
/*--------------------------------------------------------------------------*/

#include "VrySmplP.h"

#include <queue>
#include <algorithm>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace VS01P_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------ AUXILIARY FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template<typename T>
static inline T abs( const T a )
{
 return( a > 0 ? a : -a ); 
 }

/*--------------------------------------------------------------------------*/
/*------------------------- AUXILIARY CLASSES ------------------------------*/
/*--------------------------------------------------------------------------*/

namespace VS01P_di_unipi_it
{
 // this has to be defined as VS01P_di_unipi_it::ETNode instead of ::ETNode
 // for compatibility with the friend declaration in the class

 // no forward declaration is necessary due to the friend declaration
 // in the class

struct ETNode {                     // a node in the enumeration tree
 VerySimple01Problem::Weight Val;   // objective function value of the
                                    // solution corresponding to the node
 VerySimple01Problem::Index FFree;  // index (in the order by | w[ i ] | of
                                    // the first free variable, i.e., the
                                    // first (...) FFree - 1 variables are
                                    // fixed
 ETNode *Dad;                       // father in the enumeration tree
 ETNode *Nxt;                       // next node in the singly-linked list
                                    // of ETNs (for cleanup purposes)
 };
};

/*--------------------------------------------------------------------------*/

struct myLess1 {
 // comparison operator for ordering ETNodes in nonincreasing order of
 // their Val

 bool operator()( ETNode *x , ETNode *y )
 {
  return( y->Val > x->Val );
  }
 };

/*--------------------------------------------------------------------------*/

struct myLess2 {
 // comparison operator for ordering variables in nondecreasing order of
 // | w[ i ] |

 myLess2( const VerySimple01Problem::Weight *csts ) { wght = csts; }

 bool operator()( const VerySimple01Problem::Index x ,
		  const VerySimple01Problem::Index y ) const
 {
  return( abs( wght[ x ] ) < abs( wght[ y ] ) );
  }

 const VerySimple01Problem::Weight *wght;
 };

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

VerySimple01Problem::VerySimple01Problem( Index n )
{
 nvar = n;

 w = NULL;
 ord = NULL;
 Q = NULL;
 CUList = NULL;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void VerySimple01Problem::SetWeights( const Weight *wght )
{
 if( w )
  cleanup();

 w = wght;
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

void VerySimple01Problem::SolveVS01P( void )
{
 if( ! w )
  throw VS01Pexception( "SolveVS01P(): called with no weights." );

 nsol = 0;
 curr = 0;
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

VerySimple01Problem::Weight VerySimple01Problem::GetVal( void )
{
 Weight v = 0;

 if( ! nsol )  // the optimal objective function value- - - - - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < nvar ; i++ )
   if( w[ i ] > 0 )
    v += w[ i ];

  OptVal = v;
  }
 else          // another solution- - - - - - - - - - - - - - - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  priority_queue< ETNode* , vector<ETNode*> , myLess1 > *q;

  if( nsol == 1 )  // the first (possibly) nonoptimal solution- - - - - - - -
  {
   // initialize the permutation of variables

   ord = new Index[ nvar ];

   for( Index i = 0 ; i < nvar ; i++ )
    ord[ i ] = i;

   // sort variables in nondecreasing order of | w[ i ] |

   sort( ord , ord + nvar , myLess2( w ) );

   // create Q

   q = new priority_queue< ETNode* , vector<ETNode*> , myLess1 >;
   Q = q;

   // create the root node

   curr = CUList = new ETNode;

   curr->Val = OptVal;
   curr->FFree = 0;
   curr->Dad = curr->Nxt = NULL;

   // create the first son of the root node

   curr = new ETNode;

   // the first son of the root node corresponds to the optimal solution
   // with the variable i with smallest | w[ i ] | flipped (i = ord[ 0 ]);
   // its objective function value is OptVal - | w[ i ] |, since if
   // w[ i ] < 0 (w[ i ] = - | w[ i ] |) then x^*[ i ] was 0 and has to
   // be set to 1, while if w[ i ] > 0 then x^*[ i ] was 1 and has to be
   // set to 0

   curr->Val = OptVal - abs( w[ ord[ 0 ] ] );
   curr->FFree = 1;
   curr->Dad = curr->Nxt = CUList;
   CUList = curr;
   }
  else          // any other (possibly) nonoptimal solution- - - - - - - - -
  {
   // pick the best node from Q (assumed nonempty)

   q = static_cast< priority_queue< ETNode* , vector<ETNode*> ,
                    myLess1 >* >( Q );

   if( q->empty() )
    throw VS01Pexception( "VS01P::GetSol(): all solutions seen yet." );

   curr = q->top();
   q->pop();
   }

  // the return value is the objective function value of the solution- - - - -
  // in curr

  v = curr->Val;

  // now construct the first son and the brother next to the right to - - - -
  // curr and add them to Q

  Index h = curr->FFree;

  if( h < nvar )  // ... if curr actually has a first son and a brother
  {               // next to the right, i.e., not all the variables are
                  // fixed in curr

   // in the first son, all variables are as in curr except curr->FFree
   // that is flipped w.r.t. curr (and, therefore, w.r.t. the optimal
   // solution, since all variables from curr->FFree on have the same
   // value as in the optimal solution

   ETNode *fson = new ETNode;

   fson->Val = curr->Val - abs( w[ ord[ h ] ] );
   fson->FFree = h + 1;
   fson->Dad = curr;
   fson->Nxt = CUList;
   CUList = fson;

   q->push( fson );  // insert fson in Q;

   // in the brother next to the right, all variables are as in curr->Dad
   // except curr->FFree that is flipped w.r.t. curr->Dad (and, therefore,
   // w.r.t. the optimal solution; note that curr->FFree > curr->Dad->FFree,
   // i.e., curr->FFree is not fixed in curr->Dad, and therefore it has the
   // same value as in the optimal solution

   ETNode *rbrt = new ETNode;
  
   rbrt->Val = curr->Dad->Val - abs( w[ ord[ h ] ] );
   rbrt->FFree = h + 1;
   rbrt->Dad = curr->Dad;
   rbrt->Nxt = CUList;
   CUList = rbrt;

   q->push( rbrt );  // insert fson in Q;

   }           // end( if( curr has some unfixed variable )
  }            // end( else( another solution ) ) - - - - - - - - - - - - - -
               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 nsol++;
 return( v );

 }  // end( VerySimple01Problem::GetVal )

/*--------------------------------------------------------------------------*/

void VerySimple01Problem::GetSol( ZeroOne *x )
{
 if( ! ord )  // the optimal solution- - - - - - - - - - - - - - - - - - - -
 {            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < nvar ; i++ )
   if( w[ i ] > 0 )
    x[ i ] = ZeroOne( 1 );
   else
    x[ i ] = ZeroOne( 0 );
  }
 else          // another solution- - - - - - - - - - - - - - - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ! curr )
   throw VS01Pexception( "GetSol(): called before GetVal()." );

  // compute the solution climbing up from curr to the root of the
  // enumeration tree; note that the variable flipped in a node is
  // always after (in the ordering of ord[]) the variable flipped
  // in its father, if any

  Index i = nvar;
  for( ETNode *nde = curr ; ; )
  {
   // all variables from i - 1 to nde->FFree (in the ordering of ord[])
   // have the same value as in the optimal solution

   for( Index j = nde->FFree ; j < i ; )
   {
    Index h = ord[ --i ];

    if( w[ h ] > 0 )
     x[ h ] = ZeroOne( 1 );
    else
     x[ h ] = ZeroOne( 0 );
    }

   if( nde->Dad )    // nde is not the root
    nde = nde->Dad;  // climb up
   else
    break;           // done

   // the variable nde->FFree - 1 is flipped w.r.t. the optimal solution

   Index h = ord[ --i ];

   if( w[ h ] > 0 )
    x[ h ] = ZeroOne( 0 );
   else
    x[ h ] = ZeroOne( 1 );

   }  // end( for( climbing up the enumeration tree )
  }  // end( else( another solution )
 }  // end( VerySimple01Problem::GetSol )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

VerySimple01Problem::~VerySimple01Problem()
{
 cleanup();
 }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

inline void VerySimple01Problem::cleanup( void )
{
 w = NULL;

 delete[] ord;
 ord = NULL;
 
 priority_queue< ETNode* , vector<ETNode*> , myLess1 > *q =
  static_cast< priority_queue< ETNode* , vector<ETNode*> , myLess1 >* >( Q );
 Q = NULL;

 delete q;

 for( ETNode *h = CUList ; h ; )
 {
  ETNode *nxt = h->Nxt;
  delete h;
  h = nxt;
  }

 CUList = NULL;
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File VrySmplP.C -------------------------------*/
/*--------------------------------------------------------------------------*/

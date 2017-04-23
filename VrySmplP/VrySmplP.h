/*--------------------------------------------------------------------------*/
/*----------------------- File VrySmplP.h ----------------------------------*/
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
 * \date 14 - 12 - 2003
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __VrySmplP
 #define __VrySmplP      // self-identification - endif at the end of the file

/*--------------------------------------------------------------------------*/
/*--------------------------- INCLUDES -------------------------------------*/
/*--------------------------------------------------------------------------*/

#include <exception>

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace VS01P_di_unipi_it
{
 /** @namespace VS01P_di_unipi_it
     The namespace VS01_di_unipi_it is defined to hold the
     VerySimple01Problem class and all the relative stuff. It comprises the
     namespace std. */

 using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS VerySimple01Problem --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** This class solves the Very Simple Problem in 0/1 optimization (VS01P):
    given n 0/1 variables x[ i ] and n weights w[ i ], i = 0 ... n - 1, the
    problem is
    \f[
     \max \{ \sum_{i = 0}^{n - 1} w[ i ] x[ i ] : x[ i ] \in \{ 0 , 1 \}
             \hspace{1cm} , i = 0 , ... , n - 1 \}
    \f]
    The problem is so exceedingly trivial to solve (scan the variables in any
    order, assign 1 to the variable if the weight is > 0, assign 0 if the
    weight is < 0, do whatever if the weight is zero) that one may well argue
    against the need of a class for doing that. However, this class does more
    than solving the problem: it can efficiently produce all the solutions of
    the problem (the 2^n different strings of n bits) ordered by non
    increasing objective function value.

    Of course, enumerating all these solutions can never be done efficiently
    since they are exponentially many. However, we are able to generate them
    one by one, and the effort to generate each new solution is very limited:
    generating the (k+1)-th solution only requires O(n) -- the bare minimum
    necessary to write it down. Also, only O(k) memory is required to
    generate k solutions (not counting the memory required to store the
    solution themselves, though).

    This is obtained by iteratively constructing and visiting an enumeration
    tree, whose structure is somewhat different from those used e.g. in
    enumeration algorithms. The root of the tree is (one of) the optimal
    solution(s) of the problem, say x^*. Each node of the enumeration tree 
    (not only the leaf nodes, as usual in enumeration trees) contains a
    solution of the problem. The sons of the root (imagine them pictured
    from left to right) are as follows:

    - the first node is obtained by flipping variable x[ 0 ] (giving it the
      different value with respect to x^*[ 0 ]) and solving the VS01P on all
      the other variables considering that that variable is fixed, i.e.,
      keeping all the other variables as they are in x^*;

    - the second node is obtained by keeping x[ 0 ] as in x^*, flipping
      x[ 1 ] and solving the VS01P on all the other variables considering the
      first two variables as fixed, i.e., keeping all the other variables as
      they are in x^*;

    - ...

    - the n-th node is obtained by fixing all the first n - 1 variables as
      they are in x^* and flipping x[ n - 1 ].

    The tree is then constructed by recursively iterating the process on all
    the sons. In each node, the first k variables are fixed (x[ k - 1 ] is the
    one that is flipped with respect to the father node). Then, the sons of
    that node are obtained, for i = k ... n - 1, by keeping x[ 0 ] ...
    x[ i - 1 ] as in the father node, flipping x[ i ] w.r.t. the father node
    (and, therefore, w.r.t. x^*) and solving the VS01P on all the other
    variables considering the first i + 1 variables as fixed, i.e., keeping
    all the other variables as they are in x^*.

    This tree is exponential in size (of course). However, if the variables
    are properly ordered it is possible to only explicitly construct and
    visit O(k) nodes of the tree to enumerate the first k solutions (in
    nonincreasing order of the objective function value) of VS01P. The good
    order of the variables is that for nondecreasing value of | w[ i ] |;
    with this order, flipping a variable i produces "less damage" (decreases
    less the objective function value) than flipping a variable j with j > i.

    If the variables are ordered this way, then the solutions of the problem
    can be listed, in the desired order, as follows. The set of unvisited
    nodes of the tree from which the exploration must proceed, Q, is 
    initialized with the root (containing the optimal solution x^*). Each
    time the next best solution is required, the best node (the one with
    largest objective function value of the associated solution) is extracted
    from Q and the associated solution is returned. Then, the first son of
    that node is added to Q. Also, if the node has a father (it is not the
    root) it surely is the rightmost son of his father inserted in Q as yet;
    then, its brother next to the right (if any) is also added to Q. No
    other nodes need to be added to Q, since all other sons, and all other
    brothers more to the right w.r.t. the immedaite one, surely contain
    solutions with worse (not better) objective function value w.r.t. the
    two generated ones.

    This exploration strategy ensures that producing the k best solutions
    to the problem requires generating, and inserting in Q, at most 2k nodes
    of the enumeration tree. Since each node of can be represented with O(1)
    memory, the storage required to hold the representation of the fragment
    of the tree generated so far is O(k). Also, generating the (k + 1)-th
    solution requires O( lg k ) for searching the best node in Q (using e.g.
    a binary heap), O(1) for generating the corresponding two new nodes and
    O(n) to write down the solution using the information stored in the
    current fragment of the tree (climbing the tree from the node to the
    root gives the list of all and only variables than need to be flipped
    w.r.t. x^*, and those variables are even ordered). Since k is no more
    than 2^n, the overall complexity of generating a new solution is O(n).
    */

class VerySimple01Problem
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */

   typedef double Weight;

/**< Type of the weights w[ i ]. By changing this definition and recompiling
   the code works with whatever base type is chosen. It may have been set as
   a template, but it seemed overkill. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   typedef unsigned int Index;

/**< Type of the indices "i" in "w[ i ]", "x[ i ]". By changing this
   definition and recompiling the code works with whatever base type is
   chosen. It may have been set as a template, but it seemed overkill. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   typedef double ZeroOne;

/**< Type of the variables x[ i ]: must be large enough to hold both values
   0 and 1 (is there any type which can not?). By changing this definition
   and recompiling the code works with whatever base type is chosen. It may
   have been set as a template, but it seemed overkill. */

/*--------------------------------------------------------------------------*/

   class VS01Pexception : public exception {

/**< Small class for exceptions. Derives from std::exception implementing
   the virtual method what() - and since what is virtual, remember to
   always catch it by reference (catch VS01Pexception &e) if you want
   the thing to work. */

   public:
    VS01Pexception( const char *const msg = 0 ) { errmsg = msg; }

    const char* what( void ) const throw () { return( errmsg ); }

   private:
    const char *errmsg;
   };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

   VerySimple01Problem( Index n );

/**< Constructor of the class. The parameter "n" is the number of variables
   in the VS01P. */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   void SetWeights( const Weight *wght );

/**< Set a new vector of weights, effectively changing all the data of the
   VS01P. The new weights are expected to be found in the first n positions
   of the vector wght. The class has the right to retain a pointer to this
   vector and keep using it until SetWeights() is called again, so the
   caller must not change or delete the vector until it is in use by the
   class. Calling again SetWeights() destroys any existing information about
   the previous set of k-optimal solutions, restarting the generating process
   anew.

   The class does not change the vector wght (the pointer is read-only).

   This method has to be called at least once if any method in the
   following sections is to be called. */

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   void SolveVS01P( void );

/**< Solves the problem (not a big deal). This also re-initialize the
   process for generating the solutions, so the next solution obtained with
   GetSol() [see below] after a call to SolveVS01P() will always be (one of)
   the optimal one(s).

   This method has to be called at least once every time the weights change
   [see SetWeights() above] is any solution for the corresponding VS01P is
   desired. */

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution(s)
    @{ */

   Weight GetVal( void );

/**< Get the objective function value of the next best solution.

   This method can be called only after SolveVS01P() [see above]. At the
   first call it returns the optimal objective function value of VS01P. At
   the subsequent calls it starts returning the objective function values of
   all the solutions to VS01P, in nonincreasing order. If more than one
   solution have the same value of the objective function, that value will
   be returned as many times as there are solutions.

   @note This method must not be called more than 2^n times, as there are
         "only" that many different solutions. */

/*--------------------------------------------------------------------------*/

   void GetSol( ZeroOne *x );

/**< Get the next best solution.

   This method and be called only after GetVal() [see above]. It returns the
   solution having the value of the objective function returned by the latest
   call to GetVal(). The solution is written in the first n positions of the
   vector x.

   Note that the "pointer in the list of the solutions" is only "moved" by
   calls to GetVal(). That is, two subsequent calls to GetSol() with no
   calls to GetVal() in between will return the same solution. Analogously,
   calling k times GetVal() with no calls to GetSol() in between amounts
   to discarding k - 1 solutions, since only the solution corresponding to
   the last return value of GetVal() can then be retrieved. */

/*--------------------------------------------------------------------------*/

   inline Index NSol( void ) const;

/**< Returns the number of solutions enumerated so far. */

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

   inline Index Getn( void ) const;

   inline const Weight* Getw( void ) const;

/* Methods for accessing the data structures of the class. 

   Getn() returns the number of variables.

   Getw() returns (a read-only pointer to) the weights vector. */

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   ~VerySimple01Problem();

/**< Destructor of the class. */

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ FRIENDS -----------------------------------*/
/*--------------------------------------------------------------------------*/

   friend struct ETNode;

/* Node of the enumeration tree. This friend declaration serves as a
   "forward declaration" that some VS01P_di_unipi_it::ETNode exists,
   thus allowing the class to hold pointers to ETNodes. */

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void cleanup( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   Index nvar;           // number of variables
   const Weight *w;      // vector of arc weights

   Weight OptVal;        // objective function value of the optimal solution

   unsigned long nsol;   // number of solutions generated so far

   Index *ord;           // vector containing the order (permutation) of the
                         // variables (names) in nondecreasing sense of the
                         // absolute value of the corresponding w[ i ]

   void *Q;              // priority queue (heap) of the nodes of the
                         // enumeration tree created but not yet visited:
                         // it is defined "void *" to avoid having to
                         // explicitly show the "less than" function in the
                         // header file

   struct ETNode *curr;  // node corresponding to the current solution

   struct ETNode *CUList;// head of the singly-linked list containing all the
                         // nodes of the enumeration tree created so far (for
                         // cleanup purposes)

/*--------------------------------------------------------------------------*/

 };  // end( class VerySimple01Problem )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline VerySimple01Problem::Index VerySimple01Problem::NSol( void ) const
{
 return( nsol );
 }

/*--------------------------------------------------------------------------*/

inline VerySimple01Problem::Index VerySimple01Problem::Getn( void ) const
{
 return( nvar );
 }

/*--------------------------------------------------------------------------*/

inline const VerySimple01Problem::Weight* VerySimple01Problem::Getw( void )
       const
{
 return( w );
 }

/*--------------------------------------------------------------------------*/

 };  // end( namespace VS01P_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* VrySmplP.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File VrySmplP.h -------------------------------*/
/*--------------------------------------------------------------------------*/

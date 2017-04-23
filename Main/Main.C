/*--------------------------------------------------------------------------*/
/*----------------------------- File Main.C --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- main() for testing the VerySimple01Problem class.                    --*/
/*--                                                                      --*/
/*--                          VERSION 1.01                                --*/
/*--                         14 - 12 - 2003                               --*/
/*--                                                                      --*/
/*--                Original Idea and Implementation by:                  --*/
/*--                                                                      --*/
/*--                        Antonio Frangioni                             --*/
/*--                     Operations Research Group                        --*/
/*--                    Dipartimento di Informatica                       --*/
/*--                       Universita' di Pisa                            --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "VrySmplP.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <ctime>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace VS01P_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

const VerySimple01Problem::Index TOOMANY = 21;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // read command line parameters- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( argc < 3 ) {
  cerr << "Usage: " << argv[ 0 ] << " <num var> <num sol>" << endl;
  return( 1 );
  }

 VerySimple01Problem::Index nvar;
 Str2Sthg( argv[ 1 ] , nvar );

 unsigned long int nsol;
 Str2Sthg( argv[ 2 ] , nsol );

 // enter the try-block - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 try {
  // construct the VerySimple01Problem object - - - - - - - - - - - - - - - -

  VerySimple01Problem vsp( nvar );

  // create random weights in [-100, 100) - - - - - - - - - - - - - - - - - -

  VerySimple01Problem::Weight *w = new VerySimple01Problem::Weight[ nvar ];

  srand( time( 0 ) );
  for( VerySimple01Problem::Index i = 0 ; i < nvar ; i++ )
   w[ i ] = rand() % 200 - 100;

  // pass the weights - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  vsp.SetWeights( w );

  // solve the problem- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  vsp.SolveVS01P();

  // now collect nsol solutions (ordered by o.f. value) - - - - - - - - - - -

  VerySimple01Problem::ZeroOne *x = new VerySimple01Problem::ZeroOne[ nvar ];

  // print out weights

  if( nvar < TOOMANY ) {
   cout << "\tw = [" ;

   for( VerySimple01Problem::Index i = 0 ; i < nvar - 1 ; i++ )
    cout << w[ i ] << ", ";

   cout << w[ nvar - 1 ] << "]" << endl;
   }

 VerySimple01Problem::Weight pv;

 for( unsigned long int h = 0 ; h < nsol ; h++ ) {
   // get next value
   VerySimple01Problem::Weight v = vsp.GetVal();

   // get the corresponding solution
   vsp.GetSol( x );

   // check that the value is right
   VerySimple01Problem::Weight tv = 0;
   for( VerySimple01Problem::Index i = 0 ; i < nvar ; i++ )
    tv += w[ i ] * x[ i ];

   if( tv != v )
    cout << "Error: tv - v = " << tv - v << " != 0" << endl;

   // check that the order is right
   if( h )
    if( pv < v )
     cout << "Error: out-of-order solutions (previous = " << pv
	  << ", v = " << v << ")" << endl;

   pv = v;

   // print value and solution
   if( nvar < TOOMANY ) {
    cout << "v = " << v << "\tx = [";

    for( VerySimple01Problem::Index i = 0 ; i < nvar - 1 ; i++ )
     cout << x[ i ] << ", ";

    cout << x[ nvar - 1 ] << "]" << endl;
    }
   }  // end( collectiong solutions loop )

  // cleanup- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  delete[] x;
  delete[] w;

  }  // end( try-block )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // managing exceptions - - - - - - - - - - - - - - - - - - - - - - - - -
 catch( exception &e ) {
  cerr << e.what() << endl;
  return( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  return( 1 );
  }

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Main.C -------------------------------*/
/*--------------------------------------------------------------------------*/

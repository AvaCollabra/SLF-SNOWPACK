#ifndef  __SOLVER_H__
#define  __SOLVER_H__
/* **********************************************************************************************/
/*                                        stand-alone                                          */
/*                               Derived from RESEARCH VERSION 9.0                             */
/* **********************************************************************************************/
/* **********************************************************************************/
/*  Copyright WSL Institute for Snow and Avalanche Research    SLF-DAVOS           */
/* **********************************************************************************/
/* This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

/*
* GENERAL INFO
* AUTHOR: GUIDO SARTORIS  ETH ZUERICH
*/

/**
 * @file Solver.h
 * This module define the interface functions to access the direct solver and explain the
 * sequence of function calls requested to solve a linear system of equations. Solves the
 * linear equation system: [A]{X} = {B} --> {X} := [A]-1{B}, where [A] is a sparse symmetric
 * or unsymmetric, structure-symmetric NxN matrix, {B} is a single right-hand-side vector and
 * {X} is the sough solution vector both of dimension N. In the solution process the
 * right-hand-side vector {B} is overwritten by the solution vector {X}
 * When the linear system of equations stam from a FE method for a vector field or a similar
 * one is possible and advantageous to define a multiplicity factor M which is nothing but the
 * dimension of the vector field. In this case the linear system of equations is defined as we
 * would work with a scalar field instead of a vector field. In this case each defined
 * equation and unknown is replaced by a set of M equations and M unknowns and each
 * coefficients of the matrix [A] is replaced by a full (MxM) matrix. Thus, if a the
 * multiplicity factor greather than one has been defined the true number of equations and
 * unknonw is given by the product (NxM). E.g. if you want to solve a linear system of
 * dimension DIM with a full matrix [A], simply define only one equation N=1 with multiplicity
 * M=DIM
 * The algorithm is based on the LU-factorisation of the matrix [A]=[L][U] with permutations
 * of the rows and columns of [A] such as to minimize fill-in acording to the
 * multiple-minimum-degree algorithm described by: [Alan Geroge, Joseph W.H. Liu, The
 * Evolution of the minimum degree ordering algorithm, SIAM Review, Vol. 31, No. 1 pp 1-19,
 * March 1989].
 * All direct-solver (ds-)functions return 0 = FALSE if succesfull, respectively 1 = TRUE with
 * an error message on standard output if an error has occurred.
 */



/*
 * GLOBAL MATRIX DATA
 */

/*
 * START MACRO DEFINITIONS AND TYPEDEFS
 */
#define   GD_INIT(VALUE)
#define NORM_R double

#define  GD_MEM_ERR( POINTER, MSG0, MSG )                                                      \
{                                                                                              \
	if ( POINTER  ) {                                                                      \
		gd_MemErr = false;                                                             \
	} else {                                                                               \
		gd_MemErr = true; fprintf(stderr, "\n+++++ %s: %s\n", MSG0,MSG);               \
	}                                                                                      \
}                                                                                              \

#define GD_MALLOC( POINTER, TYPE, N, MSG )                                                     \
{                                                                                              \
	POINTER = (TYPE *)malloc( sizeof(TYPE)*(N+1) );                                        \
   	GD_MEM_ERR( POINTER, "NO SPACE TO ALLOCATE", MSG );                                    \
}

#define GD_REALLOC( POINTER, TYPE, N, MSG )                                                    \
{                                                                                              \
	if ( POINTER )  {                                                                      \
  		POINTER = (TYPE *)realloc( (char*)POINTER, sizeof(TYPE)*(N+1) );               \
                     GD_MEM_ERR( POINTER, "NO SPACE TO REALLOCATE", MSG );                     \
	} else {                                                                               \
  		GD_MALLOC(  POINTER, TYPE, N, MSG );                                           \
	}                                                                                      \
}

#define GD_FREE( POINTER )                                                                     \
{                                                                                              \
	if ( POINTER ) {                                                                       \
   		free ( (char*) POINTER );                                                      \
   		POINTER = NULL;                                                                \
	}                                                                                      \
}                                                                                              \

#define GD_BCOPY( FROM, TO, NBYTES )  bcopy((char*)FROM,(char*)TO,(int)(NBYTES))


/*
 * This section contains macros which are high vectorizable. On some computers they can be
 * substitued by appropriated function calls ( BLAS routines )
 */

#define VD_PRODUCT( N, X, Y, Z )                                                               \
{  double  *x_, *y_, *z_;                                                              \
   int  k_;                                                                           \
   for( x_=X, y_=Y, z_=Z, k_=N; 0<k_--; ) (*z_++) = (*x_++) * (*y_++);                         \
}

#define VD_ADD_PRODUCT( N, X, Y, Z )                                                           \
{  double  *x_, *y_, *z_;                                                              \
   int  k_;                                                                           \
   for( x_=X, y_=Y, z_=Z, k_=N; 0<k_--; ) (*z_++) += (*x_++) * (*y_++);                        \
}

#define VD_DOT(N, X, Y, RESULT )                                                               \
{  double  *x_, *y_, r_;                                                               \
   int  k_;                                                                           \
   for( x_=X, y_=Y,  r_=0.0, k_=N; 0<k_--; ) r_ += (*x_++)*(*y_++);                            \
   RESULT = r_;                                                                                \
}

#define VD_DOT3(N, X, Y, Z, RESULT )                                                           \
{  double  *x_, *y_,*z_,  r_;                                                          \
   int  k_;                                                                           \
   for( x_=X, y_=Y, z_=Z,  r_=0.0, k_=N; 0<k_--; ) r_ += (*x_++)*(*y_++)*(*z_++);              \
   RESULT = r_;                                                                                \
}

#define VD_SUM( N, X, RESULT )                                                                 \
{  double  s_, *x_;                                                                    \
   int  k_;                                                                           \
   for( s_=0.0, x_=X, k_=N; 0<k_--; ) s_ += (*x_++);                                           \
   RESULT = s_;                                                                                \
}

#define VD_SCALE( N, A, X )                                                                    \
{  double a_, *x_;                                                                     \
   int  k_;                                                                           \
   for( a_=A, x_=X, k_=N; 0<k_--; )  (*x_++) *= a_;                                            \
}

#define VD_COPY( N, A, X )                                                                     \
{  double  a_, *x_;                                                                    \
   int  k_;                                                                           \
   for( a_=A, x_=X, k_=N; 0<k_--; )  (*x_++) = a_;                                             \
}

#define VD_SCALE_SET( N, A, X, Y )                                                             \
{  double  a_, *x_, *y_;                                                               \
   int  k_;                                                                           \
   for( a_=A, x_=X, y_=Y, k_=N; 0<k_--; )  (*y_++) = a_ * (*x_++);                             \
}

#define VD_AXPY(N, A, X, Y )      /* Y[] = A*X[] + Y[] */                                      \
{  double  a_, *x_, *y_;                                                               \
   int  k_;                                                                           \
   for (x_=X, y_=Y, a_=A, k_=N; 0<k_--; ) *y_++ += (a_)*(*x_++) ;                              \
}

#define VD_AXPY_JUMP(N_B, N, JUMP, A, X, Y ) /* Y[] = A*X[] + Y[] BLOCK-WISE IN Y */           \
{  double  a_, *x_, *y_;                                                               \
   int  n_, k_;                                                                       \
   for (x_=X, y_=Y+JUMP[0], a_=A, n_=0; n_<N_B; y_+= JUMP[++n_])                               \
   for (k_=N[n_]; 0<k_--;  ) *y_++ += (a_)*(*x_++) ;                                           \
}

#define VD_AXPY_POS(N_B, N, POS, A, X, Y ) /* Y[] = A*X[] + Y[] BLOCK-WISE IN Y */             \
{  double  a_, *x_, *y_;                                                               \
   int  n_, k_;                                                                       \
   for (x_=X, y_=Y+POS[0], a_=A, n_=0; n_<N_B; y_= Y+POS[++n_])                                \
   for (k_=N[n_]; 0<k_--;  ) *y_++ += (a_)*(*x_++) ;                                           \
}

#define VD_DOT_POS(N_B, N, POS, X, Y, RESULT) /* RESULT = X[]*Y[]  BLOCK-WISE IN Y */          \
{  double  r_, *x_, *y_;                                                               \
   int  n_, k_;                                                                       \
   for (x_=X, y_=Y+POS[0], r_=0.0, n_=0; n_<N_B; y_= Y+POS[++n_])                              \
   for (k_=N[n_]; 0<k_--;  ) r_ += (*x_++)*(*y_++) ;                                           \
   RESULT = r_;                                                                                \
}

#define VD_DOT_SPARSE( N, X, INDX, Y, RESULT)                                                  \
{  int    k_, *indx_;                                                                 \
   double  r_, *x_, *y_;                                                               \
   x_ = X;  y_ = Y;  indx_ = INDX;                                                             \
   for ( r_=0.0, k_=N; 0<k_--; )  r_ += (*x_++) * y_[ (*indx_++) ];                            \
   RESULT = r_;                                                                                \
}

#define VD_AXPY_SPARSE( N, A, X, INDX, Y)    /* Y[] = A*X[] + Y[] */                           \
{  int    k_, *indx_;                                                                 \
   double  a_, *x_, *y_;                                                               \
   a_ = A;  x_ = X;  y_ = Y;  indx_ = INDX;                                                    \
   for ( k_=N; 0<k_--; )   y_[ (*indx_++) ] += a_ * (*x_++) ;                                  \
}

/*
 * END VectorData definitions
 */

/**
 * @struct SD_CHUNK_DATA
* @brief The CHUNK_DATA is used to keep track of allocated memory chunks ( block of memory ). For
* each allocated chunk we store the pointer returned from the memory allocator, so that a
* later deallocation of the memory will be possible.
*/

#define SD_CHUNK_REALLOC_SIZE  250

typedef struct
{
	int     nChunks;
	int     pChunksSize;
	char   **pChunks;
	int     TotChunkSize;

} SD_CHUNK_DATA;

#define SD_ALLOC_CHUNK(CHUNK,SIZE)                                                             \
{  if ( CHUNK.nChunks >= CHUNK.pChunksSize )                                                   \
   {  CHUNK.pChunksSize += SD_CHUNK_REALLOC_SIZE;                                              \
      GD_REALLOC( CHUNK.pChunks, char*, CHUNK.pChunksSize, "Chunk pointer Data" );             \
   }                                                                                           \
   GD_MALLOC( CHUNK.pChunks[ CHUNK.nChunks ], char, SIZE, "Chunk Data" );                      \
   CHUNK.TotChunkSize += SIZE;                                                                 \
   CHUNK.nChunks++;                                                                            \
}

#define SD_DESTROY_CHUNK(CHUNK)                                                                \
{  int i_;                                                                                     \
   for(i_=0; i_<CHUNK.nChunks; i_++) GD_FREE(CHUNK.pChunks[i_]);                               \
   GD_FREE(CHUNK.pChunks);                                                                     \
   CHUNK.TotChunkSize = 0;                                                                     \
}

/*
* BLOCK MATRIX DATA USED FOR NUMERICAL FACTORIZATION
*/

/**
 * @struct SD_ROW_BLOCK_DATA
* @brief The data structure to store the matrix for numerical factorization is a simple one. The
* matrix structure is after the mmd sorting algorithm and the symbolic factorization mainly
* composed of clustered non-zero matrix coefficients which form blocks. In this case we use a
* data strucutre to represent these row and column blocks.
* NOTE: The row blocks are simply what in the literature is specified as supernodes. We have
* kept the data structure as simple as possible to minimize the numerical operations and so
* the execution time.
* NOTE: When we define a system we have the possibility to define a multiplicity factor. This
* allows us to perform all initialization tasks ( input of incidences, symbolic
* factorization, block format calculation ) with all indices modulo the multiplicity factor
* to save time and memory. Thus the computed block format is always the same for any
* specified multiplicity factor only the block bounds are different.
* ATTENTION: The size of the permutation vector is only Dim/(Multiplicity Factor).
*/

#define  double double
#define  SD_MARKED    (1<<30) /* An flag bit used for permutation. Use:(1<<15) on PC */

typedef struct
{
	int                Row0;
	int                Row1;
	int                nCol;
	int                nColBlock;
	int                iColBlock;
	int                iFloat;

}  SD_ROW_BLOCK_DATA;

typedef struct
{
	int                 Dim;
	int                *pPerm;
	int                 nRowBlock;
	SD_ROW_BLOCK_DATA  *pRowBlock;
	int                 nColBlock;
	int                *pFirstColBlock;
	int                *pSizeColBlock;
	int                 SizeBlockJump;
	int                *pBlockJump;
	int                 SizeUpper;
	double              *pUpper;

}  SD_BLOCK_MATRIX_DATA;


#define SD_P_FIRST_COL_BLOCK(pMAT,pROW) (pMAT->pFirstColBlock + pROW->iColBlock)
#define SD_P_SIZE_COL_BLOCK(pMAT,pROW)  (pMAT->pSizeColBlock  + pROW->iColBlock)

#define SD_SEARCH_BLOCK_ROW(ROW, pROW_LOW, pROW_HIGH, pROW)                                    \
{  SD_ROW_BLOCK_DATA *low_, *high_, *mid_;                                                     \
   low_ = pROW_LOW; high_ = pROW_HIGH;                                                         \
   while( low_<=high_ )                                                                        \
   {  mid_ = low_ + ( high_ - low_ ) / 2;                                                      \
      if      ( ROW < mid_->Row0 ) high_ = mid_ - 1;                                           \
      else if ( ROW > mid_->Row1 ) low_  = mid_ + 1;                                           \
      else { pROW=mid_;  break;  }                                                             \
   }                                                                                           \
}

/*
* SPARSE MATRIX DATA USED INITIALLY TO READ THE CONNECTIVITY MATRIX
*/

// #define SPARSE_BINARY_TREE

#ifdef  DEBUG_SOLVER
	#define SD_FILLED     ( 1<<30 )
	#define SD_DEBUG_BIT  ( SD_FILLED )
	#define SD_COL(pCOL)  ( (pCOL)->Col & (~SD_DEBUG_BIT) )
#else
	#define SD_FILLED     0
	#define SD_COL(pCOL)  ( (pCOL)->Col )
#endif

typedef struct SD_COL_DATA
{
	int                  Col;

	#ifdef SPARSE_BINARY_TREE
	struct SD_COL_DATA  *Right, *Left;
	#else
	struct SD_COL_DATA  *Next;
	#endif

} SD_COL_DATA;

typedef struct SD_ROW_DATA
{
	SD_COL_DATA  *Col;

}  SD_ROW_DATA;

typedef struct
{
	int                    nRow;
	int                   *pPerm;
	int                   *pPermInv;
	int                    nSupernode;
	int                   *pSupernode;

	SD_ROW_DATA           *pRow;
	SD_CHUNK_DATA          PoolCol;
	SD_COL_DATA           *FreeCol;
	int                    nFreeCol;
	int                    nCol;

}  SD_CON_MATRIX_DATA;

#define SD_ROW_NUM(pROW, pMAT)    ( pROW - (pMAT)->pRow  )
#define SD_ROW(     NUM, pMAT)    ( (pMAT)->pRow[NUM]    )

/*
* BLOCK SPARSE MATRIX DATA USED FOR SYMBOLIC FACTORIZATION
*/

/*
* When we are building the block matrix format we need some other data structure than the
* final one. In this case the intermediary data strcuture is of smaller size than the final
* one and has a common data field. In this case with a union definition we force both data
* structures to have the same size, and the same common field so that to avoid unnecessary
* reallocation and movement of data.
*/

typedef struct SD_COL_BLOCK_DATA
{
	int Col0, Col1;
	struct SD_COL_BLOCK_DATA  *Next;

} SD_COL_BLOCK_DATA;

typedef union
{
	SD_ROW_BLOCK_DATA     UnusedData;
	struct
	{
		int                Row0;
		int                Row1;
	}  Any;
	struct
	{
		int                Row0;
		int                Row1;
		SD_COL_BLOCK_DATA *ColBlock;
	}  Data;
}  SD_TMP_ROW_BLOCK_DATA;

typedef struct
{
	int                    nRow;
	int                   *pPerm;

	int                    nRowBlock;
	SD_TMP_ROW_BLOCK_DATA *pRowBlock;
	int                    nColBlock;
	SD_CHUNK_DATA          PoolColBlock;
	SD_COL_BLOCK_DATA     *FreeColBlock;

}  SD_TMP_CON_MATRIX_DATA;

/**
 * @struct SD_MATRIX_DATA
* @brief When the user define a matrix, the software return a pointer to an opaque type i.e. a
* pointer to void as index to reference the matrix. This pointer is actually the pointer to
* the SD_MATRIX_DATA data structure. This date structure is defined as a union of differnet
* matrix data representations, and the type of data actually stored depend on the evolution
* of the algorithn.
*/

typedef enum StateType {ConMatrix, BlockConMatrix, BlockMatrix}  StateType;

typedef  struct
{
	int   nEq;
	int   nDeletedEq;
	int   Multiplicity;

	StateType State;
	/*	enum
	{  ConMatrix,
	BlockConMatrix,
	BlockMatrix
	}  State;
	*/
	union
	{  SD_CON_MATRIX_DATA      Con;
	SD_TMP_CON_MATRIX_DATA  TmpCon;
	SD_BLOCK_MATRIX_DATA    Block;
	}  Mat;

}  SD_MATRIX_DATA;
typedef  SD_MATRIX_DATA MYTYPE;


/*
* DATA FOR COLUMN AND COLUMN BLOCK ALLOCATION
*/

#define SD_ALLOC_COL(N_COL, pMAT)                                                              \
{  SD_ALLOC_CHUNK((pMAT)->PoolCol, sizeof(SD_COL_DATA)*N_COL);                                 \
   (pMAT)->FreeCol = ( SD_COL_DATA * ) (pMAT)->PoolCol.pChunks[ (pMAT)->PoolCol.nChunks-1 ];   \
   (pMAT)->nFreeCol = N_COL;                                                                   \
}

#define SD_N_ALLOC_COL  500

#define SD_GET_COL(pCOL, pMAT)                                                                 \
{  if  ( !(pMAT)->nFreeCol )  SD_ALLOC_COL(SD_N_ALLOC_COL, pMAT);                              \
   pCOL = (pMAT)->FreeCol++; (pMAT)->nFreeCol--;                                               \
}

/*
* To little speed-up memory operations for the SD_COL_BLOCK_DATA, we only allocate chunks of
* SD_COL_BLOCK_DATA and put each cell in a linked LIFO list of free SD_COL_BLOCK_DATA. When
* we need one cell of SD_COL_BLOCK_DATA we take it from the free list, by release of the data
* we put it again in the free list. A LIFO list also called a stack is important to avoid
* eccessive scattering of data in memory.
*/

#define SD_FREE_COL_BLOCK_0(pCOL_BLOCK, pMAT)                                                  \
{  (pCOL_BLOCK)->Next = (pMAT)->FreeColBlock; (pMAT)->FreeColBlock = (pCOL_BLOCK);  }

#define SD_FREE_COL_BLOCK(pCOL_BLOCK, pMAT)                                                    \
{  SD_FREE_COL_BLOCK_0(pCOL_BLOCK, pMAT); pMAT->nColBlock--;  }

#define SD_ALLOC_COL_BLOCK(N_COL_BLOCK, pMAT)                                                  \
{  SD_COL_BLOCK_DATA *pColBlock_;                                                              \
   SD_ALLOC_CHUNK((pMAT)->PoolColBlock, sizeof(SD_COL_BLOCK_DATA)*N_COL_BLOCK);                \
   pColBlock_ = ( SD_COL_BLOCK_DATA * )                                                        \
                (pMAT)->PoolColBlock.pChunks[ (pMAT)->PoolColBlock.nChunks-1 ];                \
   for(int i_=N_COL_BLOCK; 0<i_--; pColBlock_++) SD_FREE_COL_BLOCK_0(pColBlock_, pMAT);        \
}

#define SD_N_ALLOC_COL_BLOCK  500

#define SD_GET_COL_BLOCK(pCOL_BLOCK, pMAT)                                                     \
{  if  ( !(pMAT)->FreeColBlock )  SD_ALLOC_COL_BLOCK(SD_N_ALLOC_COL_BLOCK, pMAT);              \
   pCOL_BLOCK = (pMAT)->FreeColBlock; (pMAT)->FreeColBlock = (pMAT)->FreeColBlock->Next;       \
   pMAT->nColBlock++;                                                                          \
}

/*
* COLUMN DATA MANAGEMENT
*/

/**
* @brief The SD_FIND_COL macro, accept pROOT_COL as the first node of the list to start the search
* of the node with value COL. If the node is found the variable FOUND is set to TRUE and
* ppCOL will point to this node, if not found FOUND is set to FALSE and ppCOL will point to
* the entry point.
* NOTE. For 2D meshes is better to not enable SPARSE_BINARY_TREE and use a linear list of
* column coefficients instead of a binary tree of column coefficients.
*/

#ifdef SPARSE_BINARY_TREE

#define SD_FIND_COL(pROOT_COL, COL, ppCOL, FOUND)                                              \
{  SD_COL_DATA *pC_ ;                                                                          \
   FOUND   = 0;                                                                                \
   pC_     = (pROOT_COL);                                                                      \
   ppCOL = &(pROOT_COL);                                                                       \
   while ( pC_ )                                                                               \
   {  if      ( COL < SD_COL(pC_) ) {  ppCOL = &pC_->Left;  pC_ = pC_->Left;   }               \
      else if ( COL > SD_COL(pC_) ) {  ppCOL = &pC_->Right; pC_ = pC_->Right;  }               \
      else {  FOUND = 1; break;  }                                                             \
   }                                                                                           \
}

#define SD_INSERT_COL(ppCOL, pCOL, COL)                                                        \
{  pCOL->Left = pCOL->Right  = 0;                                                              \
   pCOL->Col  = COL;                                                                           \
   *ppCOL     = pCOL;                                                                          \
}

#else

#define SD_TRAVERSE_COL(pCOL, pFUNCTION, pVOID)                                                \
 {  SD_COL_DATA *pC_; for(pC_=pCOL; pC_; pC_=pC_->Next) pFUNCTION(pC_, pVOID);  }

#define SD_FIND_COL(pROOT_COL, COL, ppCOL, FOUND)                                              \
{  SD_COL_DATA *pC_ ;                                                                          \
   FOUND   = 0;                                                                                \
   pC_     = (pROOT_COL);                                                                      \
   ppCOL   = &(pROOT_COL);                                                                     \
   while ( pC_ )                                                                               \
   {  if ( COL > SD_COL(pC_)  )  { ppCOL = &pC_->Next; pC_ = pC_->Next;  }                     \
      else { if ( COL == SD_COL(pC_) ) FOUND = 1;  break;  }                                   \
   }                                                                                           \
}

#define SD_INSERT_COL(ppCOL, pCOL, COL)                                                        \
{  pCOL->Col  = COL;                                                                           \
   pCOL->Next = *ppCOL;                                                                        \
   *ppCOL     = pCOL;                                                                          \
}

#define SD_ADD_COL_1(COL, pADD_COL)                                                            \
{  SD_COL_DATA **ppC_, *pC_;                                                                   \
   int Col_, Found_;                                                                           \
   Col_ = COL;                                                                                 \
   SD_FIND_COL_LIN(pCOL, Col_, ppC_, Found_);                                                  \
   if ( !Found_ ) {  SD_GET_COL(pC_, pMat); SD_INSERT_COL(ppC_, pC_, Col_); }                  \
}

#endif


#define ERROR_SOLVER(MSG)      { printf("++++Errror:%s:%d:%s\n", __FILE__, __LINE__, MSG); return(1); }
#define USER_ERROR(MSG) { printf(ErrMsg, MSG); return 1; }
#define EXIT(MSG)  {  printf("++++Exit:%s:%d:%s\n", __FILE__, __LINE__, MSG);  exit(1);   }


#define MAX_MULT 100 /* a very big value */
#define FOR_MULT(X) for(m=0; m<Mult; m++) {X;}


typedef struct  {
	int *pC0, *pSize;
} pBLOCK;

#define BLOCK_INIT(BLOCK,pCOL0,pSIZE) { BLOCK.pC0 = pCOL0; BLOCK.pSize = pSIZE; }
#define BLOCK_NEXT(BLOCK)             ( BLOCK.pC0++,       BLOCK.pSize++ )
#define BLOCK_C0(BLOCK)                 BLOCK.pC0[0]
#define BLOCK_C1(BLOCK)                (BLOCK.pC0[0]+BLOCK.pSize[0])
#define BLOCK_SIZE(BLOCK)               BLOCK.pSize[0]

#define pC0_FIRST_COL(  pROW)  (pMatFirstColBlock + pROW->iColBlock)
#define pSIZE_FIRST_COL(pROW)  (pMatSizeColBlock  + pROW->iColBlock)

#define FIND_COL_BLOCK(pFIRST_BLK, COL, ppBLK, FOUND)                                          \
{                                                                                              \
   FOUND      = 0;                                                                             \
   SD_COL_BLOCK_DATA *pB_ = (pFIRST_BLK);                                                      \
   int Col0_  = COL+1;                                                                         \
   int Col1_  = COL-1;                                                                         \
   while ( pB_ )                                                                               \
   {  if      ( Col1_ >  pB_->Col1  ) { ppBLK = &pB_->Next; pB_ = pB_->Next; }                 \
      else if ( Col0_ >= pB_->Col0  )                                                          \
      {  FOUND = 1;                                                                            \
         if      ( Col0_==pB_->Col0 ) pB_->Col0 = COL;                                         \
         else if ( Col1_==pB_->Col1 ) pB_->Col1 = COL;                                         \
         break;                                                                                \
      }                                                                                        \
      else   break;                                                                            \
   }                                                                                           \
}

/*
* Macros to compute the triangular factorization on a symmetric matrix stored packed row-wise
* in a one dimensional array. i.e the lower matrix coefficient are not stored. This macro is
* used to invert the pivot row block if its size is greater than 1.
*/

#define FACT_SYM_MAT(MAT,N_ROW,N_COL)                                                          \
{                                                                                              \
   if ( N_ROW>1 ) {                                                                            \
   const int m_n_1=N_COL-N_ROW+1;                                                              \
   double *Mat_k=MAT;                                                                           \
   for ( int n_k=N_COL; n_k>=m_n_1; n_k-- )                                                    \
   {  double Pivot = 1./(*Mat_k);                                                               \
      double *Mat_i = Mat_k++;                                                                  \
      for ( int n_i=n_k; n_i>m_n_1; Mat_k++ )                                                  \
      {  Mat_i += n_i--;  VD_AXPY(n_i, -(*Mat_k)*Pivot, Mat_k, Mat_i);  }                      \
      Mat_k += m_n_1 - 1;                                                                      \
   }                                                                                           \
  }                                                                                            \
}

#define FACT_SYM_MAT_UNUSED(MAT,N)                                                             \
{  int   n_k, n_i;  double Pivot, *Mat_k;                                                       \
   for ( Mat_k=MAT, n_k=M; n_k>0; n_k-- )                                                      \
   {  Pivot = 1./(*Mat_k);                                                                     \
      double *Mat_i = Mat_k++;                                                                  \
      for ( n_i=n_k; n_i>1; Mat_k++ )                                                          \
      {  Mat_i += n_i--; VD_AXPY(n_i, -(*Mat_k)*Pivot, Mat_k, Mat_i);  }                       \
   }                                                                                           \
}

/*
* This macro compute the triangular factorization for a block of rows of dimension N_PIVOT
* onto another block of rows of dimension N_ROW for a block symmetric matrix stored packed
* row-wise in a one dimensional array. Schenatically we have:
*
*    N_PIVOT            N_ROW                       N[1]            N[2]
*     <--->           <-------->                 <-------->         <-->
*     X X X  -  -  -  *  *  *  *  *  *  *        *  *  *  *         *  *
*     X X X  -  -  -  *  *  *  *  *  *  *        *  *  *  *         *  *
*     X X X  -  -  -  *  *  *  *  *  *  *        *  *  *  *         *  *
*     <--------------><----------------->
*        TOT_ROW            N_COL         JUMP[1]           JUMP[2]
*                                        <------>          <------->
*                     *  *  *  *  *  *  *  -  -  *  *  *  *  -   -  *  *
*                        *  *  *  *  *  *  -  -  *  *  *  *  -   -  *  *
*                           *  *  *  *  *  -  -  *  *  *  *  -   -  *  *
*                              *  *  *  *  -  -  *  *  *  *  -   -  *  *
*
* The first values of both row block are specified by MAT0 and MAT1, N[0] and JUMP[0] are
* unused. ATTENTION: This macro change the value of N[0] which is first set to N_COL and then
* is changed continously.
*/
#if 1
#define FACT_SYM_MAT_BLOCK(N_PIVOT,TOT_ROW,N_ROW,N_COL,MAT0,DIM0,MAT1,DIM1,N_BLOCK,N,JUMP)     \
{  int n_k, i_, k__, dim_i;  double Pivot, *Mat_k0, *Mat_k, *Mat_i;                              \
   for ( Mat_k0=MAT0, k__=0, n_k=N_PIVOT; n_k>0; n_k--, k__++ )                                  \
   {  Pivot = 1./(*Mat_k0);                                                                    \
      Mat_k = Mat_k0 + TOT_ROW - k__;                                                           \
      Mat_i = MAT1;                                                                            \
      dim_i = DIM1;                                                                            \
      N[0]  = N_COL;                                                                           \
      for ( i_=N_ROW; i_>0; Mat_k++, i_-- )                                                    \
      {  VD_AXPY_JUMP(N_BLOCK, N, JUMP, -(*Mat_k)*Pivot, Mat_k, Mat_i );                       \
         N[0]-- ;                                                                              \
         Mat_i += (dim_i)--;                                                                   \
      }                                                                                        \
      Mat_k0  += DIM0 - k__;                                                                    \
   }                                                                                           \
}

#else
void FACT_SYM_MAT_BLOCK (int N_PIVOT, int TOT_ROW, int N_ROW, int N_COL, double *MAT0, int DIM0,
                         double *MAT1, int DIM1, int N_BLOCK, int *N, int *JUMP)
{
	int n_k, i_, k__, dim_i;
	double Pivot, *Mat_k0, *Mat_k, *Mat_i;
	for ( Mat_k0 = MAT0, k__ = 0, n_k = N_PIVOT; n_k > 0; n_k--, k__++ ) {
		Pivot = 1. / (*Mat_k0);
		Mat_k = Mat_k0 + TOT_ROW - k__;
		Mat_i = MAT1;
		dim_i = DIM1;
		N[0]  = N_COL;
		for ( i_ = N_ROW; i_ > 0; Mat_k++, i_-- ) {
			VD_AXPY_JUMP(N_BLOCK, N, JUMP, -(*Mat_k)*Pivot, Mat_k, Mat_i );
			N[0]-- ;
			Mat_i += (dim_i)--;
		}
		Mat_k0  += DIM0 - k__;
	}
}
#endif
/*
* This macro performs a binary search for the row: ROW. The block containing this row is
* returned by pROW
*/

#define FIRST_BLOCK_ROW(pMAT) ( pMAT->pRowBlock )
#define LAST_BLOCK_ROW(pMAT)  ( pMAT->pRowBlock + pMAT->nRowBlock - 1 )

#define SEARCH_ROW(ROW, pROW_LOW, pROW_HIGH, pROW)                                             \
{  SD_ROW_BLOCK_DATA *low_, *high_, *mid_;                                                     \
   low_ = pROW_LOW; high_ = pROW_HIGH;                                                         \
   while( low_<=high_ )                                                                        \
   {  mid_ = low_ + ( high_ - low_ ) / 2;                                                      \
      if      ( ROW < mid_->Row0 ) high_ = mid_ - 1;                                           \
      else if ( ROW > mid_->Row1 ) low_  = mid_ + 1;                                           \
      else { pROW=mid_;  break;  }                                                             \
   }                                                                                           \
}

/*
* This macro compute the block jump offsets between two rows. The first row must be a subset
* of the second one i.e. all coefficients of the first row must be present in the second one.
*/
#if 1
#define BLOCK_JUMP(nCOL0, pCOL0, pSIZE0, pCOL1, pSIZE1, JUMP)                                  \
{  int i_, *pCol0_, *pCol1_, *pSize0_,  *pSize1_, Size_, Col1_0_, Col1_1_;                     \
   pCol0_  = pCOL0;                                                                            \
   pCol1_  = pCOL1;                                                                            \
   pSize0_ = pSIZE0;                                                                           \
   pSize1_ = pSIZE1;                                                                           \
   Col1_0_ = pCol1_[0];                                                                        \
   Col1_1_ = Col1_0_ + pSize1_[0];                                                             \
   for(i_=0; i_<nCOL0; i_++)                                                                   \
   {  Size_   = 0;                                                                             \
      while( pCol1_[0] + pSize1_[0] < pCol0_[0] )                                              \
      {  Size_ += Col1_1_ - Col1_0_;                                                           \
         Col1_1_  = ( Col1_0_ = (++pCol1_)[0] ) + (++pSize1_)[0];   }                          \
      JUMP[i_] = Size_ + pCol0_[0] - Col1_0_;                                                  \
      Col1_0_  = (pCol0_++)[0] + (pSize0_++)[0];                                               \
   }                                                                                           \
}
#else
void BLOCK_JUMP(int nCOL0, int *pCOL0, int *pSIZE0, int *pCOL1, int *pSIZE1, int *JUMP)
{
	int i_, *pCol0_, *pCol1_, *pSize0_,  *pSize1_, Size_, Col1_0_, Col1_1_;
	pCol0_  = pCOL0;
	pCol1_  = pCOL1;
	pSize0_ = pSIZE0;
	pSize1_ = pSIZE1;
	Col1_0_ = pCol1_[0];
	Col1_1_ = Col1_0_ + pSize1_[0];
	for (i_=0; i_<nCOL0; i_++) {
		Size_ = 0;
		while ( pCol1_[0] + pSize1_[0] < pCol0_[0] ) {
			Size_ += Col1_1_ - Col1_0_;
			Col1_1_  = ( Col1_0_ = (++pCol1_)[0] ) + (++pSize1_)[0];
		}
		JUMP[i_] = Size_ + pCol0_[0] - Col1_0_;
		Col1_0_  = (pCol0_++)[0] + (pSize0_++)[0];
	}
}
#endif

/*
* This macro compute for a matrix stored packed row-wise in a one dimensional array the
* position of a diagonal element in a given row.
*/

#define DIAGONAL(DIM,K) ( (K)*(DIM) -( (K)*((K)-1) )/2 )

/*
* A linear search is performed in the row pROW to find the column COL. This macro use the
* column value of the next column block to determine in which column block the column is to
* be found. In this case the dimension of the search array is set to the number of column
* block minus one. If the column block is not found, the block can only be the last defined
* column block. NOTE: Here we are forced to perform a linear search because we have to
* compute the total number of column coefficients defined prior the founded column block. A
* binary search could be used if instead of the column block size we store the sum of defined
* column coefficients. This is of course possible and only little change in the software are
* necessary, however, in this case we can no more pack in one integer the data for a column
* block definition.
*/

#define SEARCH_COL(COL, ROW, pMAT, pROW, FOUND, OFFSET)                                        \
{  int *col_, *size_;                                                   \
   col_     = SD_P_FIRST_COL_BLOCK(pMAT,pROW);                                                 \
   size_    = SD_P_SIZE_COL_BLOCK( pMAT,pROW);                                                 \
   const int delta_   = ROW - pROW->Row0;                                                                \
   OFFSET   = pROW->iFloat + DIAGONAL(pROW->nCol, delta_);                                     \
   {  ++col_;                                                                                  \
      for(int i_=pROW->nColBlock-1; (i_--)>0; OFFSET += size_[0], col_++, size_++)                 \
      {  if ( COL < col_[0] )  break;  }                                                       \
       --col_;                                                                                 \
      if ( COL >= col_[0]+size_[0] ) {  FOUND = 0; }                                           \
      else                           {  FOUND = 1;  OFFSET += COL - col_[0] - delta_;  }       \
   }                                                                                           \
}



/*
 * END MACRO DEFINITIONS
 */

/**
 * @brief This is the first function to be called prior to begin to work with any matrix. The user
 * must gives the number of equations i.e. the dimension of the matrix, the type of matrix
 * i.e. a symmetric or unsymmetric matrix and the multiplicity. The multiplicity specifies
 * that each coefficient of the matrix is actually a square matrix and each equation and
 * unknown is actually a set of equations and unknowns all of dimension equal to the defined
 * multiplicity factor. Thus the true dimension of the linear system is given by the specified
 * matrix dimension times the multiplicity. Of course a multiplicity value of one is the most
 * general case, but sometines especially by vector field computations the multiplicity is
 * simply given by the dimension of the vector to be computed and this allow to speed up many
 * integer operations. Only define a multiplicity factor greather than one, if the vector
 * field components, or a subset, are fully coupled to eachother. In fact, by definition of a
 * multiplicity we always reserve memory for the full coupled system among each vector
 * component, thus the memory requirement increase quadratically with the defined multiplicity
 * factor. The function return a pointer to an opaque data type as matrix identifier.
 *
 * @param MatDim : (int) dimension of the matrix [A]. The true number of equations
 * and unknowns is given by: MatDim * Multiplicity
 * @param Multiplicity : (int) Multiplicity factor with value >= 1
 * @param **ppMat : (MYTYPE) A pointer to an opaque data type storing data related to the matrix [A]
 */
int ds_Initialize(int MatDim, int Multiplicity, MYTYPE **ppMat);

/**
 * @brief This function is needed for defining the system (matrix) connectivity i.e. the non-zero
 * coefficients of the matrix [A] i.e. which equation is connected to which one. For each
 * (finite) element we have to specifies a list of equations. Here, we assume that all
 * equations in the list are connected to eachother and thus lead to non-zero coefficients in
 * the matrix i.e. the element list of equations form a crique
 * This function is generally called for a single element (nEl = 1) or for a group of nEl
 * elements having criques of equal dimension. To contemporary define the element connectivity
 * for more elements set nEl>1 and store the list equations for each element in the array:
 * Eq[][Dim]. Of course nEl<=Dim must holds.
 * If a multiplicity factor greather than 1 has been defined we have only to define the
 * element connectivity for the first representative equations or first vector field component
 * i.e. as we would do for a scalar field, and thus all equations in the list must have a
 * number less than the specified matrix dimension
 * E.g. for an element with 4 nodes and multiplicity 3 we have a total true number of 12
 * equations forming a crique but only 4 = ( 12 / 3 ) equations must be specified and the
 * values must not exceed the specified dimension of the matrix [A].
 * After the list of equations for all elements have been specified the matrix [A] should be
 * symbolically factorized by calling the function: ds_Solve(SymbolicFactorize, ... ).
 * NOTE: Except the definition of the multiplicity in ds_Initialize(), all steps performed to
 * define the structure of matrix [A] are stricktly independent from the multiplicity
 *
 * @param *pMat0 (MYTYPE): INPUT: Pointer to the matrix [A] opaque data returned by ds_Initialize()
 * @param nEq (int) INPUT: No. of equations for one element forming a crique
 * @param Eq[] (int) INPUT: Element list of equations for more elements with equal no. of eqs.
 * @param nEl (int) INPUT: No. of elements  ( 0 <= i "<" nEq ;  0 <= e "<" nEl )
 * @param Dim (int) INPUT: first dimension of the 2D-array Eq[][Dim]
 */

int ds_DefineConnectivity( MYTYPE *pMat0, int nEq, int Eq[], int nEl, int Dim );

/**
* @brief This function assemble the element square matrix [ElMat] for one element with nEq*M x nEq*M
* real coefficients in to the global matrix [A]. If a multiplicity factor M greather than 1
* has been defined the numerical values in the element matrix [ElMat] must be forall vector
* components clustered together. E.g. for a multiplicity of 3 i.e. a 3D vector field, the 3x3
* left-upper submatrix of [ElMat] must represent the coupling between the 3 vector field
* components.
* To perform this task we also newly require the element connectivity. The list of element
* equations should be the same or a subset of the list used previously when the matrix
* connectivity has been defined with a call to: ds_DefineConnectivity().
* ATTENTION: no error is detected if some of the the element connectivity defined when
* calling ds_AssembleLocalMatrix() are not included in those previously specified when
* calling ds_DefineConnectivity()
* If the matrix [A] has been declared as symmetric only the upper triangular part of
* [ElMat], i.e. only ElMat[i][j] = ElMat[Dim*i+j] with i<=j and i,j = 0...M*nEq-1 are used
* and need to be defined. In the unsymmetric case all M*nEq x M*nEq coefficients are used. It
* is assumed that in the calling program the array [ElMat] is dimensioned as ElMat[...][Dim]
* NOTICE: Of course the parameter Dim can be greater than M*nEq: Dim >= M*nEq
* After all element matrices have been assembled the matrix [A] = [L][U] can be factorised in
* to the lower [L] and upper [U] tringular matrices by calling ds_Solve(NumericFactorize,
* ).
 * @param *pMat (MYTYPE) INPUT: pointer to the matrix [A] opaque data returned by ds_Initialize()
 * @param nEq (int) INPUT: no. of equations for one element forming a crique
 * @param Eq[] (int) INPUT: Element list of equations for one element.
 * @param Dim (int) INPUT: first dimension of the 2D-array ElMat[][Dim]
 * @param *ElMat (double) INPUT: element square matrix to be assembled in the matrix [A]
*/

int ds_AssembleMatrix(MYTYPE *pMat, int nEq, int Eq[], int Dim, double *ElMat);

/**
* The next function represents the computational kernel of this direct solver. Its
* functionality depends on the input parameter: Code whose meanings are given below.
* Generally once the matrix as been defined and the element connectivity assembled, the user
* first perform a symbolic factorization process. After assembling the numerical matrix, a
* numerical factorization step follows together with a back- and for-ward substitution to
* compute a numerical solution of the linear system. The matrix data can be resetted in order
* to solve other linear systems having the same matrix connectivity without the need to newly
* symbolic factorize the matrix.
* NOTE: If a multiplicity factor has been defined we assume that the vector components are
* clustered together in the vector pVec.
*/

typedef enum SD_MATRIX_WHAT
{
   /*
    * For symbolically factorizing and optimizing the structure of the [L] and [U] matrices
    * after the connectivity of the matrix [A] has been defined by calling
    * ds_DefineConnectivity (...)
    */
   SymbolicFactorize = 1,
   /*
    * For numerically factorizing the matrix [A] = [L][U] in to the lower resp. upper
    * triangular matrices [L] resp. [U] after [A] has been symbolically factorized by calling
    * ds_Solve(SymbolicFactorize,...) and assembled by calling ds_AssemblLocalMatrix(...) for
    * each element.
    */
   NumericFactorize  = 2,
   /*
    * For solving for a new right hand side load vector {B} after the matrix [A] has been
    * numerically factorized by calling ds_Solve(NumericFactorize,...)
    */
   BackForwardSubst  = 4,
   /*
    * For both numerically factorizing the matrix [A] and solving for a 1st right hand side
    * vector {B} i.e. is equivalent to the calls ds_Solve(NumericFactorize,...) &
    * ds_Solve(BackForwardSubst,...)
    */
   ComputeSolution   = NumericFactorize | BackForwardSubst,
   /*
    * For resetting all real coefficients of the matrix [A] to zero before reassembling a new
    * matrix [A] with identical connectivity by calling ds_AssembleLocalMatrix(...) with
    * changed local element matrices [ElMat] but the same element equations list Eq[...]
    */
   ResetMatrixData   = 8,
   /*
    * For releasing all storage space needed for solving the current problem defined by the
    * data pointed to by *pMat
    */
   ReleaseMatrixData = 16

} SD_MATRIX_WHAT;

/**
 * @param Code (SD_MATRIX_WHAT) INPUT: functionlaity code defined above
 * @param *pMat (MYTYPE) INPUT: pointer to the matrix [A] opaque data
 * @param *pX (double) INOUT: right hand side vector {B} to be overwritten by the solution vector {X}:  B[i] := X[i]
 */
int ds_Solve(SD_MATRIX_WHAT Code, MYTYPE *pMat, double *pX);

double BlockMatrixElement
  (  /* SD_BLOCK_MATRIX_DATA * pMat ,
     int Row ,
     int Col */
  );
int PrintNumMatrix
  ( /*  SD_BLOCK_MATRIX_DATA * pMat ,
     int Permuted */
  );
int TestBlockFormat
  (  /* SD_BLOCK_MATRIX_DATA * pMat */
  );
int TestNumericalFact
  (  /* SD_BLOCK_MATRIX_DATA * pMat ,
     int StoreNewData  */
  );
int ds_AssembleMatrix
  (  /* SD_MATRIX_DATA * pMat0 ,
     int nEq ,
     int Eq [ ] ,
     int Dim ,
     double * ElMat */
  );
int ds_AssembleRowCoeff
  (  SD_MATRIX_DATA * pMat0 ,
     int Row ,
     int nCol ,
     int * pCol ,
     double Diagonal ,
     double * Coeff
  );
int Permute
  (  int N ,
     int * Perm ,
     double * Vector
  );
int PermuteWithMult
  (  int N ,
     int Mult ,
     int * Perm ,
     double * Vector
  );
int InvertMatrix
  (  SD_BLOCK_MATRIX_DATA * pMat
  );
int MatrixVector
  (  SD_BLOCK_MATRIX_DATA * pMat ,
     double * X ,
     double * Y
  );
int InverseMatrixVector
  (  SD_BLOCK_MATRIX_DATA * pMat ,
     double * X
  );

//For SymFact:
int PrintInitialConFormat
  (  /* SD_CON_MATRIX_DATA * pMat ,
     char * Msg */
  );
int PrintTmpConFormat
  ( /* SD_TMP_CON_MATRIX_DATA * pMat ,
     char * Msg */
  );
int PrintConFormat
  (  /* SD_BLOCK_MATRIX_DATA * pMat ,
     char * Msg */
  );
int PSInitialConFormat
  (  /* char * FileName ,
     SD_CON_MATRIX_DATA * pMat ,
     int Permuted */
  );
int PSTmpConFormat
  (  /* char * FileName ,
     SD_TMP_CON_MATRIX_DATA * pMat */
  );
int AssembleConnectivityPermuted
  (  /* SD_CON_MATRIX_DATA * pMat ,
     int * Perm ,
     int nInc ,
     int * Inc  */
  );
int SimpleSymbolicFact
  (  /* SD_CON_MATRIX_DATA * pMat */
  );
int CheckCon
  (  /* SD_TMP_CON_MATRIX_DATA * pMat */
  );
int CheckFillIn
  (  /* SD_TMP_CON_MATRIX_DATA * pMat */
  );
int CheckEqualMat
  ( /*  SD_CON_MATRIX_DATA * pMat0 ,
     SD_TMP_CON_MATRIX_DATA * pMat1 */
  );
int CheckSupernode
  (  /* SD_CON_MATRIX_DATA * pMat ,
     int * pSize ,
     int * pPermInv */
  );
int BuildSparseConFormat
  (  /* SD_CON_MATRIX_DATA * pMat ,
     int * pRowStart0 ,
     int * pColumn0  */
  );
int ComputePermutation
  (  /* SD_CON_MATRIX_DATA * pMat */
  );
int ComputeTmpConMatrix
  (  /* SD_CON_MATRIX_DATA * pMat0 ,
     SD_TMP_CON_MATRIX_DATA * pMat  */
  );
int ComputeFillIn
  (  /* SD_TMP_CON_MATRIX_DATA * pMat */
  );
int ComputeBlockMatrix
  ( /*  SD_TMP_CON_MATRIX_DATA * pTmpMat ,
     SD_BLOCK_MATRIX_DATA * pMat ,
     int Mult */
  );
int AllocateConData
  (  int Dim ,
     SD_CON_MATRIX_DATA * pMat
  );
int ReleaseConMatrix
  (   SD_CON_MATRIX_DATA * pMat
  );
int ReleaseBlockMatrix
  ( SD_BLOCK_MATRIX_DATA * pMat
  );
int ds_DefineConnectivity
  (  /* SD_MATRIX_DATA * pMat0 ,
     int nEq ,
     int Eq [ ] ,
     int nEl ,
     int Dim  */
  );
int SymbolicFact
  ( SD_MATRIX_DATA * pMat
  );



#endif

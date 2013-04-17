/*
* GENERAL INFO
* AUTHOR: GUIDO SARTORIS ETH ZUERICH
*/

/*
* In this module are defined the function the user call
* in order to solve a linear system of equations
*/

#include <snowpack/Solver.h>

/*
* DEFINE STATEMENTS
*/
static char ErrMsg[] = "++++Errror:gs_SolveMatrix:%s\n";
static bool gd_MemErr;

/*
 * INTERFACE FUNCITONS TO ACCESS THE SOLVER
 */
int ds_Initialize(int MatDim, int Multiplicity, MYTYPE **ppMat)
{
	MYTYPE  *pMat = NULL;

	if ( Multiplicity<=0 ){
		 Multiplicity=1;
	}

	GD_MALLOC( pMat, MYTYPE, 1, "Matrix Data");
	memset( pMat,0,sizeof(MYTYPE) );
	pMat->nEq = MatDim * Multiplicity;
	pMat->Multiplicity = Multiplicity;
	if ( AllocateConData( MatDim, &pMat->Mat.Con ) )
		 return 1;

	pMat->State = ConMatrix;

	*ppMat = pMat;

	return 0;

}  /* ds_Initialize */


int ds_Solve( SD_MATRIX_WHAT Code, MYTYPE *pMat, double *X)
{
	// SymbolicFactorize

	if ( Code & SymbolicFactorize ){
		if ( Code & NumericFactorize ){
			USER_ERROR("You cannot invert the matrix symbolically and numerically contemporary");
		}

		if ( pMat->State != ConMatrix ){
			USER_ERROR("Bad Matrix Format for Symbolic Factorization");
		}

      		SymbolicFact(pMat);
	}

	// NumericFactoriz

	if ( Code & NumericFactorize ){
		if (  pMat->State != BlockMatrix ){
			USER_ERROR("Bad Matrix Format for Numerical Factorization");
		}
		InvertMatrix( &pMat->Mat.Block );
		// PrintNumMatrix(&pMat->Mat.Block,1);
	}

	// BackForwardSubst
	if ( Code & BackForwardSubst ){
		int DimTot, Mult, i;
		if (  pMat->State != BlockMatrix ){
			USER_ERROR("Bad Matrix Format for Back- For-ward Substitution");
		}
		DimTot = pMat->Mat.Block.Dim + pMat->nDeletedEq;
		Mult   = pMat->Multiplicity;

		#define PRINT_VEC(N,VEC){ int i; for(i=0; i<N; i++) printf("%e ", VEC[i]); printf("\n\n");  }

       		if ( Mult==1 ){
			Permute( DimTot, pMat->Mat.Block.pPerm, X );
		} else{
			PermuteWithMult( DimTot/Mult, Mult, pMat->Mat.Block.pPerm, X );
		}
       		InverseMatrixVector( &pMat->Mat.Block, X );
       		for(i=pMat->Mat.Block.Dim; i<DimTot; i++){
			X[i] = 0.;
		}
       		if ( Mult==1 ){
			Permute( DimTot, pMat->Mat.Block.pPerm, X );
		} else{
			PermuteWithMult( DimTot/Mult, Mult, pMat->Mat.Block.pPerm, X );
		}
	}

	// ResetMatrixData
   	if ( Code & ResetMatrixData ){
		if ( Code != ResetMatrixData ){
			USER_ERROR("You cannot reset the matrix together with other operations");
		}

      		if ( pMat->State != BlockMatrix ){
			USER_ERROR("Bad Matrix Format to reset matrix");
		}

       		memset( pMat->Mat.Block.pUpper, 0, pMat->Mat.Block.SizeUpper * sizeof(double) );
	}

   	// ReleaseMatrixData
   	if ( Code & ReleaseMatrixData ){
		if ( pMat->State == ConMatrix ){
			ReleaseConMatrix(&pMat->Mat.Con);
		} else if ( pMat->State == BlockMatrix  ){
			ReleaseBlockMatrix(&pMat->Mat.Block);
		} else ERROR_SOLVER("Unknown matrix state");{
			GD_FREE(pMat);
		}
	}

	return 0;

}  /* ds_Solve */


int ds_InitializeBoeing(int MatDim, int *pxConCon, double *pData, MYTYPE **ppMat)
{
	int  Row, Col;

	ds_Initialize(MatDim, 0, ppMat);

	for(Row=0; Row<MatDim; Row++){
		int nCol, *pCol, Inc[2];
		Inc[0] = Row;
		nCol = pxConCon[Row+1] - pxConCon[Row];
		pCol = pxConCon + pxConCon[Row] - 1 ;
		for(Col=0; Col<nCol; Col++){
			Inc[1] = pCol[Col] - 1;
			ds_DefineConnectivity( *ppMat, 2, Inc, 1, 0);
		}
	}

	SymbolicFact(*ppMat);

	for(Row=0; Row<MatDim; Row++){
		int nCol, *pCol;
		nCol = pxConCon[Row+1] - pxConCon[Row];
		pCol = pxConCon + pxConCon[Row] - 1 ;
		ds_AssembleRowCoeff(*ppMat, Row, nCol, pCol, pData[Row], pData + pxConCon[Row] - 1 );
	}

	return 0;

}  /* ds_InitializeBoeing */


int ds_MatrixConnectivity( MYTYPE *pMat0, int *pMatDim, int **ppxConCon, int *pSize)
{
	SD_CON_MATRIX_DATA   *pMat;
	int Row, Col, *pRowStart, *pColumn;
	SD_ROW_DATA *pRow;

	pMat   = &pMat0->Mat.Con;
	*pMatDim  = pMat->nRow;
	*pSize = pMat->nCol/2 + *pMatDim + 1;
	GD_MALLOC( *ppxConCon,  int, *pSize , "connectivity vector");
	if ( gd_MemErr ){
		return 1;
	}

	for (Row = 0, pRow = pMat->pRow, pRowStart = *ppxConCon, pColumn = *ppxConCon + *pMatDim + 1,
		*pRowStart = *pMatDim + 1; Row<pMat->nRow; Row++, pRow++, pRowStart++){
		int nCol;
		SD_COL_DATA *pCol;
		nCol = 0;
		pCol = pRow->Col;
		while( pCol ){
			Col =  SD_COL(pCol);
			if ( Col > Row ){
				*pColumn++ = Col;  nCol++;
			}
			pCol = pCol->Next;
		}
		pRowStart[1] = pRowStart[0] + nCol;
	}

	return 0;

}  /* ds_MatrixConnectivity */




/*
 * Beginning of NumFact.c
 */
/*
* NUMERICAL FACTORIZATION ROUTINES
*/


/**
 * @brief This function assemble the element matrix for one element and must be called for each
 * (finite) element after the element connectivity have been assembled and the matrix symbolic
 * factorized. To perform this task we also newly require the element incidences. The
 * variable: Dim specifies the dimension of the matrix: Mat which is not required to be equal
 * to the numer of element incidences: nEq.
 * ATTENTION: This function do not generate a run time error if the specified incidences have
 * not been previously defined.
 * NOTE: If the matrix has been specified as symmetric we always use only the upper part of
 * the element matrix.
 * @param *pMat0 SD_MATRIX_DATA
 * @param nEq int
 * @param Eq int
 * @param Dim int
 * @param *ElMat double
 * @return int
 */
int ds_AssembleMatrix(SD_MATRIX_DATA *pMat0, int nEq, int Eq[], int Dim, double *ElMat)
{
	SD_BLOCK_MATRIX_DATA *pMat=NULL;
	SD_ROW_BLOCK_DATA *pRow=NULL;
	int Row, Col, PermRow, PermCol, Found, Index;
	int Mult;

	pMat = &pMat0->Mat.Block;
	Mult = pMat0->Multiplicity;

	for (Row = 0; Row < nEq*Mult; Row++) {
		// PermRow = Mult*pMat->pPerm[ Eq[Row%nEq ] ] + Row/nEq;
		PermRow = Mult * pMat->pPerm[ Eq[Row/Mult] ] + Row%Mult;
		SEARCH_ROW(PermRow, FIRST_BLOCK_ROW(pMat), LAST_BLOCK_ROW(pMat), pRow);
		for (Col = 0; Col < nEq*Mult; Col++) {
			// PermCol = Mult*pMat->pPerm[ Eq[Col%nEq]  ] + Col/nEq;
			PermCol = Mult * pMat->pPerm[ Eq[Col / Mult] ] + Col%Mult;
			if ( PermCol < PermRow ) {
				continue;
			}
			SEARCH_COL(PermCol, PermRow, pMat, pRow, Found, Index);
			if ( Found ) {
				if ( Row<Col ) {
					pMat->pUpper[Index] += ElMat[ Row*Dim + Col ];
				} else {
					pMat->pUpper[Index] += ElMat[ Col*Dim + Row ];
				}
			}
		}
	}
	return 0;

}  /* ds_AssembleMatrix */


/**
 * @brief This function assemble in a given row the value of the column coefficients. ATTENTION: Is
 * supposed that the array Col use the FORTRAN notation.
 * @param *pMat0 SD_MATRIX_DATA
 * @param Row int
 * @param nCol int
 * @param *pCol int
 * @param Diagonal double
 * @param *Coeff double
 * @return int
*/
int ds_AssembleRowCoeff(SD_MATRIX_DATA * pMat0, int Row, int nCol, int * pCol, double Diagonal, double * Coeff)
{
	SD_BLOCK_MATRIX_DATA *pMat=NULL;
	SD_ROW_BLOCK_DATA *pRow=NULL;
	int Col, PermRow, PermCol, Found, Index;

	pMat = &pMat0->Mat.Block;

	PermRow =  pMat->pPerm[ Row ];
	SEARCH_ROW(PermRow, FIRST_BLOCK_ROW(pMat), LAST_BLOCK_ROW(pMat), pRow);

	SEARCH_COL(PermRow, PermRow, pMat, pRow, Found, Index);
	if ( Found ) {
		pMat->pUpper[Index] += Diagonal;
	}
	for (Col = 0; Col < nCol; Col++) {
		PermCol =  pMat->pPerm[ pCol[Col]-1 ];
		if ( PermRow < PermCol ) {
			SEARCH_COL(PermCol, PermRow, pMat, pRow, Found, Index);
			if ( Found ) {
				pMat->pUpper[Index] += Coeff[Col];
			}
		} else {
			SD_ROW_BLOCK_DATA  *pRow1=NULL;
			SEARCH_ROW(PermCol, FIRST_BLOCK_ROW(pMat), LAST_BLOCK_ROW(pMat), pRow1);
			SEARCH_COL(PermRow, PermCol, pMat, pRow1, Found, Index);
			if ( Found ) {
				pMat->pUpper[Index] += Coeff[Col];
			}
		}
	}
	return 0;

}  // ds_AssembleRowCoeff


/**
 * @brief This function permute a vector, for a given permutation vector and compute the inverse
 * permutation vector. Of course this is a trivial task if we have a second vector to store
 * the new permuted vector, but we do not want to allocated extra memory so that we use a
 * somewhat expensive algorithm which do not yet require a second storage vector. Here we
 * assume that the dimension of permutation array is less than 2**31, because we use the 31th
 * bit of each index of the permutation vector to store a flag, which tells us that we have
 * already permuted that element. The algorithm is very simple, we take one element, look if
 * we have already permuted that element, and if not we permute all elements connected to this
 * element by a permutation cycle i.e. we shift right by one all elements of that cycle.
 * ATTENTION: This function not only permute the given vector, but also compute the inverse
 * permutation vector which is stored at the place of the permutation vector.
 * @param N int
 * @param *Perm int
 * @param *Vector double
 * @return int
 */
int Permute(int N, int * Perm, double * Vector)
{
	int   i, To, From, ToNext;
	double ValueTo, Value;

	for (i = 0; i < N;  Perm[i++] &= (~SD_MARKED) ) {
		if ( Perm[i] & SD_MARKED ) {
			continue;
		}
		ValueTo = Vector[i];
		From = i;
		To = Perm[i];
		while (1) {  // follows the cycle, until we find an already permuted element
			Value      = Vector[To];
			Vector[To] = ValueTo;
			ToNext     = Perm[To];
			Perm  [To] = From | SD_MARKED;
			if ( ToNext & SD_MARKED ) {
				break;
			}
			ValueTo    = Value;
			From       = To;
			To         = ToNext;
		}
	}
	return 0;

}  // Permute


/**
 * @brief The same as before but with a multyplicity factor. We define a second function to be little
 * faster when Mult == 1. NOTE: We are limited to a maximal multiplicity factor.
 * @param N int
 * @param Mult int
 * @param *Perm int
 * @param *Vector double
 * @return int
 */
int PermuteWithMult(int N, int Mult, int *Perm, double *Vector)
{
	int   m, i, To, From, ToNext;
	double    ValueTo[MAX_MULT], Value[MAX_MULT];

	if ( Mult > MAX_MULT ) {
		printf ( "+++++ Multiplicy factor %d to large\n", Mult );
		return 1;
	}

	for (i = 0; i < N;  Perm[i++] &= (~SD_MARKED) ) {
		if ( Perm[i] & SD_MARKED )  {
			continue;
		}
		FOR_MULT( ValueTo[m] = Vector[Mult*i+m] );
		From    = i;
		To      = Perm[i];

		while (1) { // follows the cycle, until we find an already permuted element
			FOR_MULT( Value[m] = Vector[Mult*To+m] ) ;
			FOR_MULT( Vector[Mult*To+m] = ValueTo[m] );
			ToNext = Perm[To];
			Perm[To] = From | SD_MARKED;
			if ( ToNext & SD_MARKED ) {
				break;
			}
			FOR_MULT( ValueTo[m] = Value[m] );
			From = To;
			To = ToNext;
		}
	}
	return 0;

}  // PermuteWithMult


/**
 * @brief This function is the kernel of the solution algorithm, and compute the LU triangular
 * factorization of a block matrix where the blocks are so defined that there is no more
 * fill-in during the factorizatio process. The LU factorization of the matrix substitute the
 * original matrix data. The LU factorization of a matrix is very similar to the inverse of
 * the matrix and from this fact we have derived the name for this function.
 * @param *pMat SD_BLOCK_MATRIX_DATA
 * @return int
*/
int  InvertMatrix( SD_BLOCK_MATRIX_DATA *pMat )
{
	SD_ROW_BLOCK_DATA *pPivotRow, *pSearchRow, *pLastSearchRow;
	int nPivotRow,*pBlockJump, *pMatFirstColBlock, *pMatSizeColBlock;
	double *Upper;

	pMatFirstColBlock = pMat->pFirstColBlock;
	pMatSizeColBlock = pMat->pSizeColBlock  ;
	pLastSearchRow = pMat->pRowBlock + pMat->nRowBlock - 1;
	pBlockJump = pMat->pBlockJump;
	Upper = pMat->pUpper;

	/*
	* This is the pivot block loop. For each row block defined in the matrix which are
	* considered sequentially we set the selected row block to be the pivot block. If the
	* pivot block has a size which is greather than one we factorize the pivot block, note
	* this is always a dense matrix operation. Then for each other column block in the pivot
	* block we search for the corresponding row block in the list of all row blocks and
	* factorize the column-row block pair.
	*/
	for ( pPivotRow = pMat->pRowBlock, nPivotRow = pMat->nRowBlock; (nPivotRow--)>0; pPivotRow++) {
		int iColBlock, nColBlock, nCol, Row_i0, Row_i1, RowDelta, DimCol, DimCol0, TotRow, DimPivot, DimRow;
		double *PivotUpper, *RowUpper;
		pBLOCK pColBlock;

		// Initilize some data for this new pivot block
		pSearchRow  = pPivotRow;
		nColBlock   = pPivotRow->nColBlock;
		nCol        = pPivotRow->nCol;
		PivotUpper  = Upper + pPivotRow->iFloat;
		BLOCK_INIT(pColBlock, pC0_FIRST_COL(pPivotRow), pSIZE_FIRST_COL(pPivotRow) );

		// Factorize the Pivot row block, this is a dense matrix operation
		DimPivot = pPivotRow->Row1 - pPivotRow->Row0 + 1;
		FACT_SYM_MAT(PivotUpper, DimPivot, nCol);

		/*
		* Process each column block defined in the pivot row block. In this case for each
		* column block we have to found out the corresponding row block . This is an expensive
		* search operation and as been optimized. First we restrict to a minimum the search
		* interval, the lower bound is set to the last found row block and the upper bound is
		* set to the last one defined in the matrix. Before we start a the binary search for a
		* given col block we first check if the previous founded one, is still the right one.
		* For each column block is possible that we have to process more than one row blocks so
		* that this result in a double nested for loop.
		*/
		for ( iColBlock = nColBlock, TotRow    = DimPivot, Row_i0 = BLOCK_C0(pColBlock) + TotRow,
			Row_i1 = BLOCK_C1(pColBlock); iColBlock>0; BLOCK_NEXT(pColBlock),
			Row_i0 = BLOCK_C0(pColBlock), Row_i1 = BLOCK_C1(pColBlock), iColBlock-- ) {

			for ( ; Row_i0<Row_i1; TotRow += DimRow, Row_i0 += DimRow, BLOCK_SIZE(pColBlock) = DimCol0 ) {
				if ( Row_i0 > pSearchRow->Row1 ) {
					pSearchRow++;
					if ( Row_i0 > pSearchRow->Row1 ) {
						pSearchRow++;
						SEARCH_ROW(Row_i0, pSearchRow, pLastSearchRow, pSearchRow);
					}
				}

				RowDelta = Row_i0 - pSearchRow->Row0;
				DimCol   = Row_i1 - Row_i0;
				DimRow   = pSearchRow->Row1 - Row_i0 + 1;
				if ( DimCol < DimRow ) {
					DimRow = DimCol;
				}

				BLOCK_JUMP( iColBlock, pColBlock.pC0, pColBlock.pSize, pC0_FIRST_COL(pSearchRow),
					    pSIZE_FIRST_COL(pSearchRow), pBlockJump);
				pBlockJump[0] = 0 ;
				DimCol0       = BLOCK_SIZE(pColBlock); // save this value because it will change

				RowUpper = Upper + pSearchRow->iFloat + DIAGONAL(pSearchRow->nCol, RowDelta);
				FACT_SYM_MAT_BLOCK(DimPivot, TotRow, DimRow, DimCol, PivotUpper, nCol,
							RowUpper,   pSearchRow->nCol-RowDelta,
      							iColBlock,  pColBlock.pSize, pBlockJump);

				//  printf("Test Matrix I:\n"); PrintNumMatrix(pMat);
			}
		}
		/*
		{  static int i;
			printf("Test Matrix K: %d \n", i++); PrintNumMatrix(pMat);
		}
		*/
	}
	return 0;

} // InvertMatrix


/**
 * @brief Multiply the matrix with a vector. This function can be called at any time, but it makes
 * only sense before the matrix has been LU factorized.
 * @param *pMat SD_BLOCK_MATRIX_DATA
 * @param *X double
 * @param *Y double
 * @return int
 */
int  MatrixVector( SD_BLOCK_MATRIX_DATA *pMat, double *X, double *Y )
{
	SD_ROW_BLOCK_DATA *pPivotRow;
	double              *Upper;
	int                nPivotRow, *pMatFirstColBlock, *pMatSizeColBlock, nRow;

	pMatFirstColBlock = pMat->pFirstColBlock;
	pMatSizeColBlock  = pMat->pSizeColBlock;
	Upper             = pMat->pUpper;
	memset(Y, 0, sizeof(double) * (pMat->Dim) ) ;

	for( nRow = 0, pPivotRow = pMat->pRowBlock, nPivotRow = pMat->nRowBlock; (nPivotRow--)>0; pPivotRow++) {
		int         DimPivot,  nCol, dim;
		double      *PivotUpper, VectorDot;
		pBLOCK      pColBlock;

		DimPivot    = pPivotRow->Row1 - pPivotRow->Row0 + 1;
		nCol        = pPivotRow->nCol;
		PivotUpper  = Upper + pPivotRow->iFloat;
		BLOCK_INIT(pColBlock, pC0_FIRST_COL(pPivotRow), pSIZE_FIRST_COL(pPivotRow) );

		for (dim = 0; dim < DimPivot; dim++, nRow++)
		{
			VD_DOT_POS( pPivotRow->nColBlock, pColBlock.pSize, pColBlock.pC0, PivotUpper, X, VectorDot );
			Y[nRow] += VectorDot;
			pColBlock.pSize[0]--;
			pColBlock.pC0[0]++;
			PivotUpper++;
			VD_AXPY_POS(pPivotRow->nColBlock, pColBlock.pSize, pColBlock.pC0, X[nRow], PivotUpper, Y );
			PivotUpper +=  nCol - dim -1;
		}
		pColBlock.pSize[0] += dim;
		pColBlock.pC0[0]   -= dim;
	}

	return 0;

} // MatrixVector


/**
 * @brief This function perform the backward and forward substitution for a LU factorized block
 * matrix. Of course this operation is the same as multipying the inverse matrix with a vector
 * and inf fact the implementation does not much differs from a matrix vector multiplication.
 * @param *pMat SD_BLOCK_MATRIX_DATA
 * @param *X double
 * @return int
 */
int  InverseMatrixVector( SD_BLOCK_MATRIX_DATA *pMat, double *X )
{
	SD_ROW_BLOCK_DATA *pPivotRow;
	double *Upper;
	int nPivotRow, *pMatFirstColBlock, *pMatSizeColBlock, nRow;

	pMatFirstColBlock = pMat->pFirstColBlock;
	pMatSizeColBlock  = pMat->pSizeColBlock;
	Upper             = pMat->pUpper;

	// Forward substitution, solve the lower triangular system, where in the diagonal we have all ones.
	for( nRow = 0, pPivotRow = pMat->pRowBlock, nPivotRow = pMat->nRowBlock; (nPivotRow--)>0; pPivotRow++) {
		int         DimPivot, nCol, dim;
		double      *PivotUpper, Alpha;
		pBLOCK      pColBlock;

		DimPivot    = pPivotRow->Row1 - pPivotRow->Row0 + 1;
		nCol        = pPivotRow->nCol;
		PivotUpper  = Upper + pPivotRow->iFloat;
		BLOCK_INIT(pColBlock, pC0_FIRST_COL(pPivotRow), pSIZE_FIRST_COL(pPivotRow) );

		for (dim = 0; dim < DimPivot; dim++) {
			Alpha = -X[nRow++] / (*PivotUpper++);
			pColBlock.pSize[0]--;
			pColBlock.pC0[0]++;

			//printf("Alpha %f Mat0 %f Vec0 %f\n", Alpha, PivotUpper[0], X[pColBlock.pC0[0]]);
			VD_AXPY_POS(pPivotRow->nColBlock, pColBlock.pSize, pColBlock.pC0, Alpha, PivotUpper, X );
			PivotUpper +=  nCol - dim -1;

			// { int i; for(i=0; i<pMat->Dim; i++) printf("%f\n", X[i]); printf("\n");  }
		}
		pColBlock.pSize[0] += dim;
		pColBlock.pC0[0]   -= dim;
	}

	// Backward substitution, solve the upper triangular system
	for (nRow = pMat->Dim - 1, nPivotRow = pMat->nRowBlock,
		pPivotRow = pMat->pRowBlock + nPivotRow - 1; (nPivotRow--)>0; pPivotRow--) {
		int         DimPivot, nCol, dim;
		double      *PivotUpper, VectorDot;
		pBLOCK      pColBlock;

		DimPivot    = pPivotRow->Row1 - pPivotRow->Row0 + 1;
		nCol        = pPivotRow->nCol;
		PivotUpper  = Upper + pPivotRow->iFloat + DIAGONAL(nCol,DimPivot-1) + 1;
		BLOCK_INIT(pColBlock, pC0_FIRST_COL(pPivotRow), pSIZE_FIRST_COL(pPivotRow) );

		pColBlock.pC0[0]   += DimPivot;
		pColBlock.pSize[0] -= DimPivot;
		for(dim = DimPivot-1; dim >= 0; nRow--, dim--) {
			VD_DOT_POS(pPivotRow->nColBlock, pColBlock.pSize, pColBlock.pC0, PivotUpper, X, VectorDot );

			// printf("Dot %f Mat0 %f Vec0 %f\n", VectorDot, PivotUpper[0], X[pColBlock.pC0[0]]);

			pColBlock.pC0[0]--;
			pColBlock.pSize[0]++;
			X[nRow] -= VectorDot;
			PivotUpper--;
			X[nRow] /= PivotUpper[0];
			PivotUpper -=  nCol - dim ;

			//   { int i; for(i=0; i<pMat->Dim; i++) printf("%f\n", X[i]); printf("\n");  }
		}
	}

	return 0;

} // InverseMatrixVector

/*
 * End of NumFact.c
 */

/*
 * Beginning of SymbFact.c
 */

/*
 * MMD FUNCTIONS
 */

/**
* @brief MmdUpdate ---- multiple minimum degree update
* purpose -- this routine updates the degrees of nodes after a multiple elimination step.
* input parameters --
 * @param ehead -- int the beginning of the list of eliminated nodes (i.e., newly formed elements).
 * @param neqns -- int number of equations.
 * @param (xadj,adjncy) --int  adjacency structure.
 * @param delta -- int tolerance value for multiple elimination.
 * @param maxint -- int maximum machine representable (short) integer.
 * updated parameters --
 * @param mdeg -- int new minimum degree after degree update.
 * @param (head,forward,backward) -- int degree doubly linked structure.
 * @param qsize -- int size of supernode.
 * @param list,marker vector for degree update.
 * @param *tag -- int tag value.
*/
static void MmdUpdate( int ehead,  int neqns,  int *xadj,  int *adjncy,  int delta,  int *mdeg,
  			int *head, int *forward,  int *backward, int *qsize, int *list,  int *marker,  int maxint,  int *tag)
{
      int  deg, deg0, element, enode, fnode, i, iq2, istop,
           istart, j, jstop, jstart, link, mdeg0, mtag, nabor,
           node, q2head, qxhead;

      mdeg0 = *mdeg + delta;
      element = ehead;

n100:
	if ( element <= 0 ) {
		return;
	}

	// for each of the newly formed element, do the following.
	// reset tag value if necessary.
	mtag = *tag + mdeg0;
	if ( mtag >= maxint ) {
		*tag = 1;
		for ( i = 1; i <= neqns; i++ ) {
			if ( marker[i] < maxint ) {
				marker[i] = 0;
			}
		}
		mtag = *tag + mdeg0;
	}

	/*
	* create two linked lists from nodes associated with 'element':
	* one with two nabors (q2head) in the adjacency structure, and the
	* other with more than two nabors (qxhead). also compute 'deg0',
	* number of nodes in this element.
	*/
	q2head = 0;
	qxhead = 0;
	deg0   = 0;
	link   = element;

n400:
	istart = xadj[link];
	istop = xadj[link+1] - 1;
	for ( i = istart; i <= istop; i++ ) {
		enode = adjncy[i];
		link = -enode;
		if ( enode < 0 )  {
			goto n400;
		}
		if ( enode == 0 ) {
			break;
		}
		if ( qsize[enode] != 0 ) {
			deg0 += qsize[enode];
			marker[enode] = mtag;

			//'enode' requires a degree update
			if ( backward[enode] == 0 ) {
				// place either in qxhead or q2head list.
				if ( forward[enode] != 2 ) {
					list[enode] = qxhead;
					qxhead = enode;
				} else {
					list[enode] = q2head;
					q2head = enode;
				}
			}
		}
	}

	// for each node in q2 list, do the following.
	enode = q2head;
	iq2   = 1;

n900:
	if ( enode <= 0 ) {
		goto n1500;
	}
	if ( backward[enode] != 0 ) {
		goto n2200;
	}
	(*tag)++;
	deg = deg0;

	// identify the other adjacent element nabor.
	istart = xadj[enode];
	nabor = adjncy[istart];
	if ( nabor == element ) {
		nabor = adjncy[istart+1];
	}
	link = nabor;
	if ( forward[nabor] >= 0 ) {
		// nabor is uneliminated, increase degree count.
		deg += qsize[nabor];
		goto n2100;
	}

	// the nabor is eliminated. for each node in the 2nd element
	// do the following.
n1000:
	istart = xadj[link];
	istop = xadj[link+1] - 1;
	for ( i = istart; i <= istop; i++ ) {
		node = adjncy[i];
		link = -node;
		if ( node != enode ) {
			if ( node < 0  ) {
				goto n1000;
			}
			if ( node == 0 ) {
				goto n2100;
			}
			if ( qsize[node] != 0 ) {
				if ( marker[node] < *tag ) {
					// 'node' is not yet considered.
					marker[node] = *tag;
					deg += qsize[node];
				} else {
					if ( backward[node] == 0 ) {
						if ( forward[node] == 2 ) {
							// 'node' is indistinguishable from 'enode'.
							// merge them into a new supernode.
							qsize[enode] += qsize[node];
							qsize[node] = 0;
							marker[node] = maxint;
							forward[node] = -enode;
							backward[node] = -maxint;
						} else {
							// 'node' is outmacthed by 'enode'
							if (backward[node]==0) {
								backward[node] = -maxint;
							}
						}
					}
				}
			}
		}
	}
	goto n2100;

n1500:
	// for each 'enode' in the 'qx' list, do the following.
	enode = qxhead;
	iq2 = 0;

n1600:
	if ( enode <= 0 )  {
		goto n2300;
	}
	if ( backward[enode] != 0 ) {
		goto n2200;
	}
	(*tag)++;
	deg = deg0;

	// for each unmarked nabor of 'enode', do the following.
	istart = xadj[enode];
	istop = xadj[enode+1] - 1;
	for ( i = istart; i <= istop; i++ ) {
		nabor = adjncy[i];
		if ( nabor == 0 ) {
			break;
		}
		if ( marker[nabor] < *tag ) {
			marker[nabor] = *tag;
			link = nabor;
			if ( forward[nabor] >= 0 ) {
				// if uneliminated, include it in deg count.
				deg += qsize[nabor];
			} else {
n1700:
				// if eliminated, include unmarked nodes in this
				// element into the degree count.
				jstart = xadj[link];
				jstop = xadj[link+1] - 1;
				for ( j = jstart; j <= jstop; j++ ) {
					node = adjncy[j];
					link = -node;
					if ( node < 0 ) {
						goto n1700;
					}
					if ( node == 0 ) {
						break;
					}
					if ( marker[node] < *tag ) {
						marker[node] = *tag;
						deg += qsize[node];
					}
				}
			}
		}
	}

n2100:
	// update external degree of 'enode' in degree structure,
	// and '*mdeg' if necessary.
	deg = deg - qsize[enode] + 1;
	fnode = head[deg];
	forward[enode] = fnode;
	backward[enode] = -deg;
	if ( fnode > 0 ) {
		backward[fnode] = enode;
	}
	head[deg] = enode;
	if ( deg < *mdeg ) {
		*mdeg = deg;
	}

n2200:
	// get next enode in current element.
	enode = list[enode];
	if ( iq2 == 1 ) {
		goto n900;
	}
	goto n1600;

n2300:
	// get next element in the list.
	*tag = mtag;
	element = list[element];
	goto n100;

}  // MmdUpdate


/**
* @brief MmdElimin -- multiple minimum degree elimination
* Purpose -- This routine eliminates the node mdeg_node of minimum degree from the adjacency
* structure, which is stored in the quotient graph format. It also transforms the quotient
* graph representation of the elimination graph.
* Input parameters --
 * @param mdeg_node -- node of minimum degree.
 * @param maxint -- estimate of maximum representable (short) integer.
 * @param tag -- tag value.
* Updated parameters --
 * @param (xadj,adjncy) -- updated adjacency structure.
 * @param (head,forward,backward) -- degree doubly linked structure.
 * @param qsize -- size of supernode.
 * @param marker -- marker vector.
 * @param list -- temporary linked list of eliminated nabors.
*/
static void MmdElimin(int mdeg_node, int *xadj, int *adjncy, int *head, int *forward, int *backward,
				  int *qsize, int *list, int *marker, int maxint, int tag)
{
	int   element, i,   istop, istart, j,
		jstop, jstart, link,
		nabor, node, npv, nqnbrs, nxnode,
		pvnode, rlmt, rloc, rnode, xqnbr;

	// find the reachable set of 'mdeg_node' and
	// place it in the data structure.
	marker[mdeg_node] = tag;
	istart = xadj[mdeg_node];
	istop = xadj[mdeg_node+1] - 1;

	// 'element' points to the beginning of the list of
	// eliminated nabors of 'mdeg_node', and 'rloc' gives the
	// storage location for the next reachable node.
	element = 0;
	rloc = istart;
	rlmt = istop;
	for ( i = istart; i <= istop; i++ ) {
		nabor = adjncy[i];
		if ( nabor == 0 ) {
			break;
		}
		if ( marker[nabor] < tag ) {
			marker[nabor] = tag;
			if ( forward[nabor] < 0 ) {
				list[nabor] = element;
				element = nabor;
			} else {
				adjncy[rloc] = nabor;
				rloc++;
			}
		}
	}

	// merge with reachable nodes from generalized elements.
	while ( element > 0 ) {
		adjncy[rlmt] = -element;
		link = element;

		n400:
		jstart = xadj[link];
		jstop = xadj[link+1] - 1;
		for ( j = jstart; j <= jstop; j++ ) {
			node = adjncy[j];
			link = -node;
			if ( node < 0 ) {
				goto n400;
			}
			if ( node == 0 ) {
				break;
			}
			if ((marker[node]<tag)&&(forward[node]>=0)) {
				marker[node] = tag;
				//use storage from eliminated nodes if necessary.
				while ( rloc >= rlmt ) {
					link = -adjncy[rlmt];
					rloc = xadj[link];
					rlmt = xadj[link+1] - 1;
				}
				adjncy[rloc] = node;
				rloc++;
			}
		}
		element = list[element];
	}
	if ( rloc <= rlmt ) {
		adjncy[rloc] = 0;
	}
	// for each node in the reachable set, do the following.
	link = mdeg_node;

n1100:
	istart = xadj[link];
	istop = xadj[link+1] - 1;
	for ( i = istart; i <= istop; i++ ) {
		rnode = adjncy[i];
		link = -rnode;
		if ( rnode < 0 ) {
			goto n1100;
		}
		if ( rnode == 0 ) {
			return;
		}

		// 'rnode' is in the degree list structure.
		pvnode = backward[rnode];
		if (( pvnode != 0 ) && ( pvnode != (-maxint) )) {
			// then remove 'rnode' from the structure.
			nxnode = forward[rnode];
			if ( nxnode > 0 ) {
				backward[nxnode] = pvnode;
			}
			if ( pvnode > 0 ) {
				forward[pvnode] = nxnode;
			}
			npv = -pvnode;
			if ( pvnode < 0 ) {
				head[npv] = nxnode;
			}
		}

		// purge inactive quotient nabors of 'rnode'.
		jstart = xadj[rnode];
		jstop = xadj[rnode+1] - 1;
		xqnbr = jstart;
		for ( j = jstart; j <= jstop; j++ ) {
			nabor = adjncy[j];
			if ( nabor == 0 ) {
				break;
			}
			if ( marker[nabor] < tag ) {
				adjncy[xqnbr] = nabor;
				xqnbr++;
			}
		}

		// no active nabor after the purging.
		nqnbrs = xqnbr - jstart;
		if ( nqnbrs <= 0 ) {
			// merge 'rnode' with 'mdeg_node'.
			qsize[mdeg_node] += qsize[rnode];
			qsize[rnode] = 0;
			marker[rnode] = maxint;
			forward[rnode] = -mdeg_node;
			backward[rnode] = -maxint;
		} else {
			// flag 'rnode' for degree update, and
			// add 'mdeg_node' as a nabor of 'rnode'.
			forward[rnode] = nqnbrs + 1;
			backward[rnode] = 0;
			adjncy[xqnbr] = mdeg_node;
			xqnbr++;
			if ( xqnbr <= jstop ) {
				adjncy[xqnbr] = 0;
			}
		}
	}
	return;

} // MmdElimin

/**
* @brief MmdInit -- mult minimum degree initialization
* purpose -- this routine performs initialization for the multiple elimination version of the
*   minimum degree algorithm.
* input parameters --
 * @param neqns -- number of equations.
 * @param xadj -- adjacency structure.
* output parameters --
 * @param (head,dfrow,backward) -- degree doubly linked structure.
 * @param qsize -- size of supernode ( initialized to one).
 * @param list -- linked list.
 * @param marker -- marker vector.
*/
static int  MmdInit(int neqns, int *xadj, int *head, int *forward, int *backward, int *qsize, int *list, int *marker)
{
	int  fnode, ndeg, node;

	for ( node = 1; node <= neqns; node++ ) {
		head[node]   = 0;
		qsize[node]  = 1;
		marker[node] = 0;
		list[node]   = 0;
	}

	// initialize the degree doubly linked lists.
	for ( node = 1; node <= neqns; node++ ) {
		ndeg = xadj[node+1] - xadj[node] + 1;
		fnode = head[ndeg];
		forward[node] = fnode;
		head[ndeg] = node;
		if ( fnode > 0 ) {
			backward[fnode] = node;
		}
		backward[node] = -ndeg;
	}
	return 0;

} // MmdInit

/**
* @brief MmdNumbering -- multi minimum degree numbering
* purpose -- this routine performs the final step in producing the permutation and inverse
* permutation vectors in the multiple elimination version of the minimum degree ordering
* algorithm.
* input parameters --
 * @param neqns -- number of equations.
 * @param qsize -- size of supernodes at elimination.
* updated parameters --
 * @param invp -- inverse permutation vector. on input, if qsize[node] = 0, then node has been
 * @param merged into the node -invp[node]; otherwise,
 * @param -invp[node] is its inverse labelling.
 * output parameters --
 * @param perm -- the permutation vector.
*/
static void  MmdNumbering(int neqns, int *perm, int *invp, int *qsize, int *nsize)
{
	int father, nextf, node, nqsize, num, root;

	for ( *nsize=0, node = 1; node <= neqns; node++ ) {
		nqsize = qsize[node];
		if ( nqsize <= 0 ) {
			perm[node] = invp[node];
		} else {
			perm[node] = -invp[node];
			(*nsize)++;
		}
	}

	// for each node which has been merged, do the following.
	for ( node = 1; node <= neqns; node++ ) {
		if ( perm[node] <= 0 ) {
			// trace the merged tree until one which has not
			// been merged, call it root.
			father = node;
			while ( perm[father] <= 0 ) {
				father = - perm[father];
			}

			// number node after root.
			root = father;
			num = perm[root] + 1;
			invp[node] = -num;
			perm[root] = num;

			// shorten the merged tree.
			father = node;
			nextf = - perm[father];
			while ( nextf > 0 ) {
				perm[father] = -root;
				father = nextf;
				nextf = -perm[father];
			}
		}
	}

	// ready to compute perm. USE C NOTATION !!! MODIFIED GS
	for ( node = 1; node <= neqns; node++ ) {
		num        = -invp[node];
		invp[node] = num  - 1;
		perm[num ] = node - 1;
	}
	return;

}  // MmdNumbering

/**
* @brief RunMmd -- multiple minimum external degree
* purpose -- this routine implements the minimum degree algorithm. it makes use of the
* implicit representation of elimination graphs by quotient graphs, and the notion of
* indistinguishable nodes. It also implements the modifications by multiple elimination and
* minimum external degree. Caution -- the adjacency vector adjncy will be destroyed.
* Input parameters --
 *  @param neqns -- number of equations.
 *  @param (xadj, adjncy) -- the adjacency structure.
 *  @param delta -- tolerance value for multiple elimination.
 *  @param maxint -- maximum machine representable (short) integer (any smaller estimate will do)
*   for marking nodes.
* Output parameters --
 *  @param nsize -- number of supernodes.
 *  @param perm -- the minimum degree ordering.
 *  @param invp -- the inverse of perm.
 *  @param *ncsub -- an upper bound on the number of nonzero subscripts for the compressed storage scheme.
*  Working parameters --
 *  @param head -- vector for head of degree lists.
 *  @param invp -- used temporarily for degree forward link.
 *  @param perm -- used temporarily for degree backward link.
 *  @param qsize -- vector for size of supernodes.
 *  @param list -- vector for temporary linked lists.
 *  @param marker -- a temporary marker vector.
* Subroutines used -- MmdElimin, MmdInit, MmdNumbering, MmdUpdate.
*/
static void  RunMmd(int neqns, int *xadj, int *adjncy, int *invp, int *perm, int delta, int *head, int *qsize, int *nsize, int *list, int *marker, int maxint, int *ncsub)
{
	int  ehead, i, mdeg, mdlmt, mdeg_node, nextmd, num, tag;
	if ( neqns <= 0 ) {
		return;
	}

	// initialization for the minimum degree algorithm.
	*ncsub = 0;
	MmdInit( neqns, xadj, head, invp, perm, qsize, list, marker );

	//  'num' counts the number of ordered nodes plus 1.
	num = 1;

	// eliminate all isolated nodes.
	nextmd = head[1];
	while  ( nextmd > 0 ) {
		mdeg_node = nextmd;
		nextmd = invp[mdeg_node];
		marker[mdeg_node] = maxint;
		invp[mdeg_node] = -num;
		num = num + 1;
	}

	// search for node of the minimum degree. 'mdeg' is the current
	// minimum degree; 'tag' is used to facilitate marking nodes.
	if ( num > neqns ) {
		goto n1000;
	}
	tag = 1;
	head[1] = 0;
	mdeg = 2;

	// infinite loop here !
	while ( 1 ) {
		while ( head[mdeg] <= 0 ) {
			mdeg++;
		}

		// use value of 'delta' to set up 'mdlmt', which governs
		// when a degree update is to be performed.
		mdlmt = mdeg + delta;
		ehead = 0;

n500:
		mdeg_node = head[mdeg];
		while ( mdeg_node <= 0 ) {
			mdeg++;
			if ( mdeg > mdlmt ) {
				goto n900;
			}
			mdeg_node = head[mdeg];
		}

		//  remove 'mdeg_node' from the degree structure.
		nextmd = invp[mdeg_node];
		head[mdeg] = nextmd;
		if ( nextmd > 0 ) {
			perm[nextmd] = -mdeg;
		}
		invp[mdeg_node] = -num;
		*ncsub += mdeg + qsize[mdeg_node] - 2;
		if ( (num+qsize[mdeg_node]) > neqns ) {
			goto n1000;
		}

		//  eliminate 'mdeg_node' and perform quotient graph
		//  transformation. reset 'tag' value if necessary.
		tag++;
		if ( tag >= maxint ) {
			tag = 1;
			for ( i = 1; i <= neqns; i++ ) {
				if ( marker[i] < maxint ) {
					marker[i] = 0;
				}
			}
		}
		MmdElimin( mdeg_node, xadj, adjncy, head, invp, perm, qsize, list, marker, maxint, tag );
		num += qsize[mdeg_node];
		list[mdeg_node] = ehead;
		ehead = mdeg_node;
		if ( delta >= 0 ) {
			goto n500;
		}

n900:
		// update degrees of the nodes involved in the
		// minimum degree nodes elimination.
		if ( num > neqns ) {
			goto n1000;
		}
		MmdUpdate( ehead, neqns, xadj, adjncy, delta, &mdeg, head, invp, perm, qsize, list, marker, maxint, &tag );

	}

n1000:

	MmdNumbering ( neqns, perm, invp, qsize, nsize );

} // RunMmd

/**
* @brief This function build a compressed sparse matrix format which is the common used one. The
* value pRowStart[Row] specifies the offset in the array pColumn[] where the column
* coefficients for the row: Row are starting to be defined. The column coefficients for a row
* are always sorted. Of course the array pRowStart has a dimension of plus one with respect
* to the matrix dimension. ATTENTION: Here we use a FORTRAN notation, all values in the
* arrays pRowStart and pColumn are incremented by one.
 * @param pMat SD_CON_MATRIX_DATA
 * @param *pRowStart0 int
 * @param *pColumn0 int
 * @return int
*/
int BuildSparseConFormat(SD_CON_MATRIX_DATA *pMat, int *pRowStart0, int *pColumn0)
{
	 int          i, *pRowStart, *pColumn;
	SD_ROW_DATA *pRow;


	for (i = pMat->nRow, pRow = pMat->pRow, pRowStart = pRowStart0,
	      pColumn = pColumn0, *pRowStart =1; i>0; i--,pRow++, pRowStart++) {
		int nCol;
		SD_COL_DATA *pCol;
		nCol = 0;
		pCol = pRow->Col;
		while ( pCol ) {
			*pColumn++ = SD_COL(pCol)+1;
			pCol = pCol->Next; nCol++;
		}
		pRowStart[1] = pRowStart[0] + nCol;
	}


	return 0;

}  // BuildSparseConFormat

/**
* @brief This function compute the permutation array and thus run the mmd algorithm.
 * @param *pMat SD_CON_MATRIX_DATA
 * @return int
*/
int ComputePermutation( SD_CON_MATRIX_DATA *pMat )
{
	 int *head;     /* array 0..maxN */
	int *list;     /* array 0..maxN */
	int *marker;   /* array 0..maxN */
	int *xadj;     /* array 0..maxN */
	int *adjncy;   /* array 0..maxCol */
	int maxint = 32000;   /* use a better value */
	int ncsub;
	int delta;
	int nEq;

	// Allocate memory requested by the mmd algorithm
	nEq = pMat->nRow;

	GD_MALLOC( pMat->pSupernode, int, nEq , "supernode size vector");
	GD_MALLOC( pMat->pPerm,      int, nEq , "permutation vector");
	GD_MALLOC( pMat->pPermInv,   int, nEq , "inverse permutation vector");
	GD_MALLOC( head,             int, nEq*4 + 1 + pMat->nCol, "mmd tmp memory");
	list   = head + 1*nEq;
	marker = head + 2*nEq;
	xadj   = head + 3*nEq;
	adjncy = head + 4*nEq + 1;

	/*
	* Run the mmd algorithm. The mmd implementation is a translation of a FORTRAN version so
	* that we pass all array pointers decremented by one.
	*/
	delta = 0;   // Is this a good value?

	BuildSparseConFormat(pMat, xadj, adjncy);

	RunMmd(pMat->nRow, xadj-1, adjncy-1, pMat->pPerm-1, pMat->pPermInv-1, delta,
		head-1, pMat->pSupernode-1, &pMat->nSupernode, list-1, marker-1, maxint, &ncsub );

	// Release temporary data used only by mmd
	GD_FREE(head);


	return 0;

}  // ComputePermutation


/**
* @brief This routine set up a temporary adjacency block matrix format for the purpose to speed-up
* the symbolic factorization process and is called after the permutation vector has been
* computed. From this step all indices for rows and columns are permuted indices. The row
* blocks are a result of the mmd algorithm and so need not to be computed.
* Here we have to compute the column blocks for the upper matrix including the diagonal. Of
* course the matrix is still sparse because at this time the fill-in is still unknown but it
* is however possible that blocks with size larger than one will appear.
* NOTE: All data allocated in pMat0 is released.
 * @param *pMat0 SD_CON_MATRIX_DATA
 * @param *pMat SD_TMP_CON_MATRIX_DATA
 * @return int
*/
int ComputeTmpConMatrix(SD_CON_MATRIX_DATA *pMat0, SD_TMP_CON_MATRIX_DATA *pMat)
{
	SD_TMP_ROW_BLOCK_DATA *pRowBlock;
	SD_ROW_DATA           *pRow;
	int          PermRow, *pPerm, *pPermInv, *pSupernode;

	if ( sizeof(SD_ROW_BLOCK_DATA)     != sizeof(SD_TMP_ROW_BLOCK_DATA ) ||
		sizeof(SD_TMP_ROW_BLOCK_DATA) >  sizeof(SD_ROW_BLOCK_DATA)         ) {
		EXIT("DATA STRUCTURE INCOMPATIBILITY");
	}

	// First allocate definetively the matrix row block data
	memset( pMat, 0, sizeof(SD_TMP_CON_MATRIX_DATA) );
	pMat->nRowBlock = pMat0->nSupernode;
	GD_MALLOC( pRowBlock, SD_TMP_ROW_BLOCK_DATA, pMat->nRowBlock, "Row Block Allocation");
	pMat->nRow      = pMat0->nRow;
	pMat->pPerm     = pMat0->pPerm;
	pMat->pRowBlock = pRowBlock;

	pSupernode      = pMat0->pSupernode;
	pPerm           = pMat0->pPerm;
	pPermInv        = pMat0->pPermInv;
	pRow            = pMat0->pRow;

	/*
	* Set the temporary adjacency block data. For each row block process each single column
	* coefficients and the diagonal one.
	*/
	for (PermRow = 0; PermRow < pMat->nRow; PermRow++) {
		int Supernode, Row, PermCol, Found;
		SD_COL_DATA       *pCol;
		SD_COL_BLOCK_DATA **ppColBlock, *pColBlock, *pFreeColBlock, *pStartColBlock;

		Row                  = pPermInv  [PermRow];
		Supernode            = pSupernode[Row];
		if ( Supernode <= 0 ) {
			continue;
		}
		pRowBlock->Data.Row0  = PermRow;
		pRowBlock->Data.Row1  = PermRow + Supernode - 1;

		// Insert the diagonal block as first block
		SD_GET_COL_BLOCK(pColBlock, pMat);
		pColBlock->Next = 0;
		pColBlock->Col0 =  pColBlock->Col1 = PermRow;
		pRowBlock->Data.ColBlock = pColBlock;

		/*
		* Insert the other column blocks. Foreach defined column coefficient first discard all
		* column coefficients lower the diagonal. Then check if the last defined column block
		* can be used as first one to start the linear search.
		* ATTENTION: It is possible that if the block is found and the upper bound is
		* incremented by one that the updated block will be adjacent to the next column blocks.
		* This condition is checked and if acknolegded one block (the second one) is removed.
		* NOTE: We have not to correctly initilize the value of ppColBlock to the true first
		* pointer of the list of clumn block because if a block is not found this will never be
		* the first one.
		*/
		for ( pCol = pRow[Row].Col; pCol; pCol = pCol->Next ) {
			PermCol = pPerm[ SD_COL(pCol) ];
			if ( PermRow > PermCol ) {
				continue;
			}
			if ( PermCol < pColBlock->Col0  ) {
				pStartColBlock = pRowBlock->Data.ColBlock;
			} else {
				pStartColBlock = pColBlock;
			}
			ppColBlock = &pStartColBlock;
			FIND_COL_BLOCK(pStartColBlock, PermCol, ppColBlock, Found);
			if ( !Found ) {
				SD_GET_COL_BLOCK(pColBlock, pMat);
				pColBlock->Next = *ppColBlock;
				*ppColBlock     = pColBlock;
				pColBlock->Col0 = pColBlock->Col1 = PermCol;
			} else {
				pColBlock     = (*ppColBlock);
				pFreeColBlock = pColBlock->Next;
				if ( pFreeColBlock && ( pColBlock->Col1 +1 ) >= pFreeColBlock->Col0 ) {
					pColBlock->Col1 = pFreeColBlock->Col1;
					pColBlock->Next = pFreeColBlock->Next;
					SD_FREE_COL_BLOCK(pFreeColBlock, pMat);
				}
			}
		}

		// Process next row block
		pRowBlock++;
	}

	// Release all the previously used adjaceny matrix data
	GD_FREE(pMat0->pRow);
	GD_FREE(pMat0->pPermInv);
	GD_FREE(pMat0->pSupernode);
	SD_DESTROY_CHUNK(pMat0->PoolCol);

	return 0;

}  // ComputeTmpConMatrix


void MERGE_COL_BLOCK(SD_COL_BLOCK_DATA *pCOL0, SD_COL_BLOCK_DATA **ppCOL1, SD_TMP_CON_MATRIX_DATA *pMAT)
{
	SD_COL_BLOCK_DATA  *pUp_, *pLo_, **ppLo_, *pC_;

	for ( ppLo_ =  ppCOL1, pLo_  = *ppLo_, pUp_ = pCOL0;  pUp_ ; pUp_=pUp_->Next ) {
		// skeep all lower blocks not intersecting and ahead of the upper block
		while ( pLo_ && ( pLo_->Col1 + 1 < pUp_->Col0 ) ){
				ppLo_= &pLo_->Next;
				pLo_ = pLo_->Next;
		}

		/*
		* check for a nil lower block or a lower block non intersecting and beyond the upper one.
		* If acknowledged insert a new lower block equal to the upper one and continue with the
		* next upper block
		*/
		if ( !pLo_ || pLo_->Col0 > pUp_->Col1 + 1 ) {
			SD_GET_COL_BLOCK(pC_, pMAT);
			pC_->Col0 = pUp_->Col0 ;
			pC_->Col1 = pUp_->Col1 ;
			pC_->Next = pLo_;
			*ppLo_    = pC_;
			ppLo_     = &pC_->Next;
			continue;
		}

		// the lower and upper block are now intersecting, take the larger bounds for the lower block
		if ( pLo_->Col0 > pUp_->Col0 ) {
			pLo_->Col0 = pUp_->Col0;
		}
		if ( pLo_->Col1 < pUp_->Col1 ) {
			pLo_->Col1 = pUp_->Col1;
		}

		// free all next lower blocks if fully included in the upper one
		while ( pC_ = pLo_->Next, pC_ && ( pC_->Col1 <= pUp_->Col1 ) ) {
			pLo_->Next = pC_->Next;
			SD_FREE_COL_BLOCK(pC_, pMAT);
		}
		// free the next lower block if partially included or adjacent to the upper one
		if ( pC_ = pLo_->Next, pC_ && ( pC_->Col0 - 1 <= pUp_->Col1 ) ) {
			pLo_->Next = pC_->Next;
			pLo_->Col1 = pC_->Col1;
			SD_FREE_COL_BLOCK(pC_, pMAT);
		}
	}
}

/**
* @brief This function compute the fill-in by factorizing symbolically the matrix.
 * @param *pMat SD_TMP_CON_MATRIX_DATA
 * @return int
*/
int  ComputeFillIn(SD_TMP_CON_MATRIX_DATA *pMat)
{
	SD_ROW_BLOCK_DATA *pPivotRow, *pSearchRow, *pLastSearchRow;
	int               nPivotRow;

	#define pFIRST_COL_BLOCK(pROW_BLOCK)  ( ( (SD_TMP_ROW_BLOCK_DATA *)pROW_BLOCK)->Data.ColBlock )

	/*
	* This is the pivot block loop. For each row block defined in the matrix which are
	* considered sequentially we set the selected row block to be the pivot block.
	*/
	pPivotRow      = ( SD_ROW_BLOCK_DATA *) pMat->pRowBlock;
	pLastSearchRow = pPivotRow + pMat->nRowBlock - 1;

	for ( nPivotRow = pMat->nRowBlock; (nPivotRow--)>0; pPivotRow++) {
		int                Row_i0, Row_i1, DimRow, Col0;
		SD_COL_BLOCK_DATA  *pColBlock;

		/*
		* Process each column block defined in the pivot row block.In this case for each column
		* block we have to found out the corresponding row block. This is an expensive search
		* operation and has been optimized. First we restrict to a minimum the search interval,
		* the lower bound is set to the last found row block and the upper bound is set to the
		* last one defined in the matrix. Before we start the binary search for a given column
		* block we first check if the previous founded one, is still the right one. For each
		* column block is possible that we have to process more than one row blocks so that
		* this result in a double nested for loop.
		*/
		pSearchRow = pPivotRow;
		pColBlock  = pFIRST_COL_BLOCK(pPivotRow); /* the first pColBlock is never nil !!! */
		Row_i0     = pColBlock->Col0;
		Row_i1     = pColBlock->Col1;
		DimRow     = pPivotRow->Row1 - pPivotRow->Row0 + 1;
		goto  SkeepFirstMergeOperation;

		do {
			Row_i0  = pColBlock->Col0;
			Row_i1  = pColBlock->Col1;
			for ( ; Row_i0<=Row_i1; Row_i0 += DimRow ) {
				if ( Row_i0 > pSearchRow->Row1 ) {
					pSearchRow++;
					if ( Row_i0 > pSearchRow->Row1 ) {
						pSearchRow++;
						SD_SEARCH_BLOCK_ROW(Row_i0, pSearchRow, pLastSearchRow, pSearchRow);
					}
				}
				DimRow = pSearchRow->Row1 - Row_i0 + 1;

				Col0 = pFIRST_COL_BLOCK(pSearchRow)->Col0;
				MERGE_COL_BLOCK(pColBlock, &pFIRST_COL_BLOCK(pSearchRow), pMat);
				pFIRST_COL_BLOCK(pSearchRow)->Col0 = Col0;

				SkeepFirstMergeOperation:;
			}
		} while ( pColBlock = pColBlock->Next, pColBlock );
	}

	return 0;

} // ComputeFillIn

/**
* @brief Here, we define the matrix as used by the numerical factorization algorithm. If the
* multiplicity factor is not 1 at this time we define the true row and column blocks.
 * @param *pTmpMat SD_TMP_CON_MATRIX_DATA
 * @param *pMat SD_BLOCK_MATRIX_DATA
 * @param Mult int
 * @return int
*/
int ComputeBlockMatrix( SD_TMP_CON_MATRIX_DATA *pTmpMat, SD_BLOCK_MATRIX_DATA *pMat, int Mult)
{
	SD_TMP_ROW_BLOCK_DATA *pTmpRowBlock;
	SD_ROW_BLOCK_DATA     *pRowBlock;
	SD_COL_BLOCK_DATA     *pColBlock;
	int                    k, nTotCol, nTotColBlock, MaxColBlock, *pFirstColBlock, *pSizeColBlock;

	memset( pMat, 0, sizeof(SD_BLOCK_MATRIX_DATA) );
	pMat->Dim        = pTmpMat->nRow*Mult;
	pMat->pPerm      = pTmpMat->pPerm;
	pMat->nRowBlock  = pTmpMat->nRowBlock;
	pMat->nColBlock  = pTmpMat->nColBlock;
	pMat->pRowBlock  = (SD_ROW_BLOCK_DATA*)pTmpMat->pRowBlock;

	/*
	* Alloc memory to store new column block data. Traverse all the temporary column blocks
	* and store the data newly in compact form.
	*/
	GD_MALLOC( pMat->pFirstColBlock, int, pMat->nColBlock, "Matrix First Col. Block" );
	GD_MALLOC( pMat->pSizeColBlock,  int, pMat->nColBlock, "Matrix Size  Col. Block" );
	if ( gd_MemErr )
		return 1;

	for( k = pMat->nRowBlock, pFirstColBlock = pMat->pFirstColBlock, pSizeColBlock = pMat->pSizeColBlock, nTotCol = 0,
		nTotColBlock = 0, MaxColBlock = 0, pTmpRowBlock = pTmpMat->pRowBlock; k>0; k--, pTmpRowBlock++ )
	{
		int Delta, nCol, nColBlock;

		for ( pColBlock = pTmpRowBlock->Data.ColBlock, nCol = 0, nColBlock = 0; pColBlock;
				pColBlock = pColBlock->Next, nColBlock++, pFirstColBlock++, pSizeColBlock++ ) {
			pFirstColBlock[0] = Mult*pColBlock->Col0;
			pSizeColBlock [0] = Mult*( pColBlock->Col1 - pColBlock->Col0 + 1 );
			nCol             += pSizeColBlock[0];
		}

		pRowBlock = (SD_ROW_BLOCK_DATA*)pTmpRowBlock;
		pRowBlock->nCol      = nCol;
		pRowBlock->nColBlock = nColBlock;
		pRowBlock->iColBlock = nTotColBlock;
		pRowBlock->iFloat    = nTotCol;
		if ( MaxColBlock < nColBlock ) {
			MaxColBlock = nColBlock;
		}
		Delta            = Mult * ( pRowBlock->Row1 - pRowBlock->Row0 + 1 );
		pRowBlock->Row0 *= Mult;
		pRowBlock->Row1  = ( pRowBlock->Row1 + 1 ) * Mult - 1;
		nTotCol         += Delta * nCol - ( Delta*(Delta-1) )/2;
		nTotColBlock    += nColBlock;
	}
	if ( nTotColBlock != pTmpMat->nColBlock ) {
		EXIT("Wrong column block count");
	}

	/*
	* Release old column block data in the temporary connectivity matrix: pTmpMat but not the
	* row block data which is the same one used by the new block matrix: pMat.
	*/
	SD_DESTROY_CHUNK(pTmpMat->PoolColBlock);

	// Alloc memory to store numerical matrix data
	pMat->SizeUpper = nTotCol;
	GD_MALLOC( pMat->pUpper, double, pMat->SizeUpper, "Matrix Data" );
	memset( pMat->pUpper, 0, sizeof(double)*pMat->SizeUpper);

	pMat->SizeBlockJump = MaxColBlock;
	GD_MALLOC( pMat->pBlockJump, int, pMat->SizeBlockJump, "Jump Data" );

	if ( gd_MemErr )
		return 1;

	return 0;

}  // ComputeBlockMatrix

/*
* MATRIX DEFINITON AND ELEMENT INCIDENCES ASSEMBLING FUNCTIONS
*/

/**
* @brief Allocate the matrix row data in order to store the connectivity matrix data.
 * @param Dim int
 * @param *pMat SD_CON_MATRIX_DATA
 * @return int
*/
int AllocateConData( int Dim, SD_CON_MATRIX_DATA *pMat )
{
	pMat->nRow = Dim;
	GD_MALLOC( pMat->pRow, SD_ROW_DATA, pMat->nRow, "Row Allocation");
	memset( pMat->pRow, 0, sizeof(SD_ROW_DATA)*pMat->nRow );
	if ( gd_MemErr ) {
		ERROR_SOLVER("Memory Error");
	}

	return 0;

}  // AllocateConData

int ReleaseConMatrix( SD_CON_MATRIX_DATA *pMat )
{
	GD_FREE(pMat->pPerm);
	GD_FREE(pMat->pSupernode);
	GD_FREE(pMat->pPermInv);
	GD_FREE(pMat->pRow);
	SD_DESTROY_CHUNK(pMat->PoolCol);

	return 0;

}  // ReleaseConMatrix

/**
* @brief Release all the data allocated for the numerical factorization algorithm.
 * @param *pMat SD_BLOCK_MATRIX_DATA
 * @return int
*/
int ReleaseBlockMatrix( SD_BLOCK_MATRIX_DATA *pMat )
{
	GD_FREE(pMat->pFirstColBlock);
	GD_FREE(pMat->pSizeColBlock);
	GD_FREE(pMat->pRowBlock);
	GD_FREE(pMat->pBlockJump);
	GD_FREE(pMat->pUpper);
	GD_FREE(pMat->pPerm);

	return 0;

}  // ReleaseBlockMatrix

/**
* @brief This function assemble the element connnectivity for one or more elements in order to build
* a sparse matrix format. Of course we only store the upper part of the connectivity matrix
* because we only consider structure symmetric matrices.
 * @param *pMat0 SD_MATRIX_DATA
 * @param nEq int
 * @param Eq (int [])
 * @param nEl int
 * @param Dim int
 * @return int
*/
int ds_DefineConnectivity(SD_MATRIX_DATA *pMat0, int nEq, int Eq[], int nEl, int Dim )
{
	int e, i, j, Row_i;
	SD_ROW_DATA        *pRow_i;
	SD_CON_MATRIX_DATA *pMat;

	pMat = &pMat0->Mat.Con;

	for (e = 0; e < nEl; Eq += Dim, e++) {
		for (i = 0; i < nEq; i++) {
			Row_i  = Eq[i];
			pRow_i = &SD_ROW(Row_i, pMat);

			for (j = 0; j < nEq; j++) {
				int Col_j, Found;
				SD_COL_DATA **ppC, *pCol;
				Col_j  = Eq[j];
				if ( Row_i == Col_j ) {
					continue;
				}
				SD_FIND_COL(pRow_i->Col, Col_j, ppC, Found);
				if ( !Found ) {
					SD_GET_COL(pCol, pMat);
					SD_INSERT_COL(ppC, pCol, Col_j);
					pMat->nCol++;
					// printf("Inserting %d %d\n", Row_i, Col_j);
				}
			}
		}
	}
	return 0;

}  // ds_DefineConnectivity

int SymbolicFact(SD_MATRIX_DATA *pMat)
{
	SD_TMP_CON_MATRIX_DATA  TmpConMat;
	SD_BLOCK_MATRIX_DATA    BlockMat;

	ComputePermutation( &pMat->Mat.Con);
	ComputeTmpConMatrix(&pMat->Mat.Con, &TmpConMat);
	ComputeFillIn(&TmpConMat);
	ComputeBlockMatrix(&TmpConMat, &BlockMat, pMat->Multiplicity);
	pMat->State     = BlockMatrix;
	pMat->Mat.Block = BlockMat;

	return 0;

}  // SymbolicFact

/*
 * End of SymbFact.c
 */


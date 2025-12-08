/*
 * RMQRMM64.cpp
 *
 *  Created on: 16-06-2014
 *      Author: hector
 */

#include "includes/RMQRMM64.h"

bool RMQRMM64::TRACE = false;
bool RMQRMM64::RUNTEST = false;
bool RMQRMM64::SHOW_SIZE = false;
uint RMQRMM64::TEST = 100000;

RMQRMM64::RMQRMM64(char *fileName){
	this->BBLK = ceilingLog64(BLK, 2);
	this->N8BLK = BLK/8;
	this->MBLK = BLK/2;
	this->BMBLK = ceilingLog64(MBLK, 2);
	this->BSS = ceilingLog64(SS, 2);

	loadDS(fileName);

	if (RUNTEST){
		test_search_min_block();
		test_sumAtPos();
		test_rank_1();
		test_select_1();
		test_positionMinblock();
		if (leaves>0) test_rmqi();
	}
}

void RMQRMM64::init(ulong len) {
	this->BBLK = ceilingLog64(BLK, 2);
	this->N8BLK = BLK/8;
	this->MBLK = BLK/2;
	this->BMBLK = ceilingLog64(MBLK, 2);
	this->BSS = ceilingLog64(SS, 2);

	if (TRACE){
		cerr << "*** RMQ structure for A[0.." << len-1 << "] and leaf length S = " << BLK << endl;
		cerr << "Create BP sequence for 2D Min Heap..." << endl;
	}

	sizeRMM = 3072; // = 2*512 + 2048;  size for T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[]
	if (TRACE || SHOW_SIZE) cerr << " ** size of fixed small tables 3072 Bytes" << endl;

	nP = (len+1)<<1;
	ulong lenP = nP >> BW64;
	if (nP % W64)
		lenP++;
	P = new ulong[lenP];
	sizeRMM += lenP*sizeof(ulong);
	//cerr << "nP " << nP << endl;
	if (TRACE || SHOW_SIZE) cerr << " ** size of topology " << lenP*sizeof(ulong) << " Bytes" << endl;
}

// bitsPC are the number of bits for each cell in A[0..len-1]
// if deleteA=true then the array A will the delete after to create the BP sequence P
RMQRMM64::RMQRMM64(ulong *A, uint bitsPC, ulong len, bool deleteA){
	init(len);
	ulong i, n_i, n_top, pos;
	StackBP *Q = new StackBP();
	StackBP *R;
	long int top = -1;

	// [1] insert root in Q and put an opening parenthesis
	i = pos = 0;
	setBit64(P, pos);
	pos++;
	Q->val = top;
	Q->next = NULL;

	// insert A[0]
	R = new StackBP();
	R->val = 0;
	R->next = Q;
	Q = R;
	setBit64(P, pos);
	pos++;
	top = 0;

	for(i=1; i<len; i++){
		// [2] append a closing parenthesis for each value j > 0 stored in Q, such that A[j] ≥ A[i], deleting these values from Q.
		n_top = getNum64(A, top*bitsPC, bitsPC);
		n_i = getNum64(A, i*bitsPC, bitsPC);
		while(top>-1 && n_top >= n_i){
			cleanBit64(P, pos);
			pos++;
			R = Q;
			Q = Q->next;
			delete R;
			top = Q->val;
			if(top>-1)
				n_top = getNum64(A, top*bitsPC, bitsPC);
		}
		// [3] insert A[i] into Q and add an opening parenthesis
		R = new StackBP();
		R->val = i;
		R->next = Q;
		Q = R;
		setBit64(P, pos);
		pos++;
		top = i;
	}
	// [4] add an opening parenthesis for each value stored in Q
	while(Q){
		cleanBit64(P, pos);
		pos++;
		R = Q;
		Q = Q->next;
		delete R;
	}

	if(deleteA){
		cerr << " Deleting array A[] ... " << endl;
		delete [] A;
		cerr << " Array deleted ! " << endl;
	}

	if(pos != 2*(len+1)){
		cerr << " ERROR. parentheses created = " << pos << " != " << 2*(len+1) << endl;
		exit(0);
	}

	if (TRACE){
		ulong lenP = nP >> BW64;
		if (nP % W64)
			lenP++;
		cerr << " BP sequence, P[0.." << nP-1 << "]" << endl;
		for (ulong pos=0; pos<lenP; pos++)
			printBitsUlong(P[pos]);
		cerr << endl;
	}
	createMinMaxTree();
}

RMQRMM64::RMQRMM64(short int *A, ulong len) {
	init(len);

	ulong i, pos=0;
	StackBP *Q = new StackBP();
	StackBP *R;
	long int top = -1;

	// [1] insert root in Q and put an opening parenthesis
	setBit64(P, pos);
	pos++;
	Q->val = top;
	Q->next = NULL;
	for(i=0; i<len; i++){
		// [2] append a closing parenthesis for each value j > 0 stored in Q, such that A[j] ≥ A[i], deleting these values from Q.
		while(top>-1 && A[top] >= A[i]){
			cleanBit64(P, pos);
			pos++;
			R = Q;
			Q = Q->next;
			delete R;
			top = Q->val;
		}
		// [3] insert A[i] into Q and add an opening parenthesis
		R = new StackBP();
		R->val = i;
		R->next = Q;
		Q = R;
		setBit64(P, pos);
		pos++;
		top = i;
	}
	// [4] add an opening parenthesis for each value stored in Q
	while(Q){
		cleanBit64(P, pos);
		pos++;
		R = Q;
		Q = Q->next;
		delete R;
	}

	if(pos != 2*(len+1)){
		cerr << " ERROR. parentheses created = " << pos << " != " << 2*(len+1) << endl;
		exit(0);
	}
	if (TRACE){
		ulong lenP = nP >> BW64;
		if (nP % W64)
			lenP++;
		cerr << " BP sequence, P[0.." << nP-1 << "]" << endl;
		for (ulong pos=0; pos<lenP; pos++)
			printBitsUlong(P[pos]);
		cerr << endl;
	}
	createMinMaxTree();
}

RMQRMM64::RMQRMM64(int *A, ulong len) {
	init(len);

	ulong i, pos=0;
	StackBP *Q = new StackBP();
	StackBP *R;
	long int top = -1;

	// [1] insert root in Q and put an opening parenthesis
	setBit64(P, pos);
	pos++;
	Q->val = top;
	Q->next = NULL;
	for(i=0; i<len; i++){
		// [2] append a closing parenthesis for each value j > 0 stored in Q, such that A[j] ≥ A[i], deleting these values from Q.
		while(top>-1 && A[top] >= A[i]){
			cleanBit64(P, pos);
			pos++;
			R = Q;
			Q = Q->next;
			delete R;
			top = Q->val;
		}
		// [3] insert A[i] into Q and add an opening parenthesis
		R = new StackBP();
		R->val = i;
		R->next = Q;
		Q = R;
		setBit64(P, pos);
		pos++;
		top = i;
	}
	// [4] add an opening parenthesis for each value stored in Q
	while(Q){
		cleanBit64(P, pos);
		pos++;
		R = Q;
		Q = Q->next;
		delete R;
	}

	if(pos != 2*(len+1)){
		cerr << " ERROR. parentheses created = " << pos << " != " << 2*(len+1) << endl;
		exit(0);
	}
	if (TRACE){
		ulong lenP = nP >> BW64;
		if (nP % W64)
			lenP++;
		cerr << " BP sequence, P[0.." << nP-1 << "]" << endl;
		for (ulong pos=0; pos<lenP; pos++)
			printBitsUlong(P[pos]);
		cerr << endl;
	}
	createMinMaxTree();
}

RMQRMM64::RMQRMM64(long int *A, ulong len) {
	init(len);

	ulong i, pos=0;
	StackBP *Q = new StackBP();
	StackBP *R;
	long int top = -1;

	// [1] insert root in Q and put an opening parenthesis
	setBit64(P, pos);
	pos++;
	Q->val = top;
	Q->next = NULL;
	for(i=0; i<len; i++){
		// [2] append a closing parenthesis for each value j > 0 stored in Q, such that A[j] ≥ A[i], deleting these values from Q.
		while(top>-1 && A[top] >= A[i]){
			cleanBit64(P, pos);
			pos++;
			R = Q;
			Q = Q->next;
			delete R;
			top = Q->val;
		}
		// [3] insert A[i] into Q and add an opening parenthesis
		R = new StackBP();
		R->val = i;
		R->next = Q;
		Q = R;
		setBit64(P, pos);
		pos++;
		top = i;
	}
	// [4] add an opening parenthesis for each value stored in Q
	while(Q){
		cleanBit64(P, pos);
		pos++;
		R = Q;
		Q = Q->next;
		delete R;
	}

	if(pos != 2*(len+1)){
		cerr << " ERROR. parentheses created = " << pos << " != " << 2*(len+1) << endl;
		exit(0);
	}
	if (TRACE && false){
		ulong lenP = nP >> BW64;
		if (nP % W64)
			lenP++;
		cerr << " BP sequence, P[0.." << nP-1 << "]" << endl;
		for (ulong pos=0; pos<lenP; pos++)
			printBitsUlong(P[pos]);
		cerr << endl;
	}
	createMinMaxTree();
}

void RMQRMM64::createMinMaxTree(){
	ulong groups, leavesUp, sizeDS;
	ulong i, j, position, father, node, cont, child, rb;
	int miniBck;
	ulong *auxL, *auxR;
	int *Aux_BkM;

	if(RUNTEST){ // test of balanced sequence
		int sum = 0;
		ulong r=0;

		for (; r<nP; r++){
			if(readBit64(P, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cerr << " ERROR. P[] is not a balanced sequence of parentheses !! " << endl;
			exit(0);
		}else
			cerr << " P[] is a well balanced sequence of parentheses !! " << endl;
	}

	ulong nb = nP / BLK;
	if (nb%2 && nb > 1)
		nb--;
	nBin = nb*BLK;
	nW = nP/W64;
	if (nP%W64)
		nW++;

	if (nBin >= nP){
		nBin = nP;
		lenLast = 0;
	}else
		lenLast = nP - nBin;

	if (TRACE)
		cerr << "Create RangeMinMaxTree_Bin with length N = 2n = " << nP << ", nBin: " << nBin << ", nW: " << nW << endl;

	this->leaves = nBin/BLK;
	uint h = 0;

	if (TRACE  && false){
		ulong i;
		cerr << "___________________________________________________________" << endl;
		cerr << "P_bin :";
		for (i=0; i<nBin; i++){
			cerr << readBit64(P, i);
			if ((i+1)%BLK == 0)
				cerr << "-";
		}
		cerr << endl << "Last Block :" << endl;
		for (; i<nP; i++)
			cerr << readBit64(P, i);
		cerr << endl;
	}

	if (leaves){
		h = 1 + (uint)(log(leaves)/log(2));
		groups = (ulong)pow(2, (double)h-1.0);	// number of nodes at level (h-1)
		firstLeaf = (ulong)pow(2, (double)h) - 1;	// the left-most leaf
		leavesBottom = (leaves - groups)<<1;
		leavesUp = leaves - leavesBottom;
		cantIN = firstLeaf - leavesUp;
		cantN = cantIN + leaves;
	}else
		cantN = cantIN =firstLeaf = leavesBottom = leavesUp = 0;

	if (TRACE && leaves>0){
		cerr << "leaves: " << leaves << ", leaves Up: " << leavesUp << ", leaves Bottom: " << leavesBottom << endl;
		cerr << "internal nodes: " << cantIN << ", first leaf: " << firstLeaf << ", total nodes: " << cantN << endl;
	}

	// Create arrays of excess in blocks...
	//cerr << "Creating arrays of excess in blocks..." << endl;
	createTables();

	// Create min - max relatives values to each internal node.
	//cerr << "Createting min - max relatives values to each internal node..." << endl;

	if (leaves){
		//cerr << "Creating auxL[] and auxR[]... " << 2*cantN*sizeof(ulong) << " Bytes" << endl;
		auxL = new ulong[cantN];	// these are auxiliary vectors to facility the compute of intervals in each internal nodes. These will be delete later.
		auxR = new ulong[cantN];

		// Step 1: process the n bits in groups of size S, and create the leaves in each array:	intervals of relatives excess per block, in MinLeaf and MaxLeaf,
		//         and the relative sum of excess per block, in ELeaf
		// 'i' is for bits in P, 'cont' is for count the size of each block, and 'node' is the current leaf position to set min-max values...
		//cerr << "Step 1: process the n bits in groups of size s, and create the leaves in each array..." << endl;
		j = firstLeaf; // this is for auxiliary vectors
		ulong bitBottom = leavesBottom*BLK;
		for (i=node=cont=0; i<nBin; i++){
			if (i == bitBottom){
				node = leavesBottom;	// we move up one level
				j = cantIN;
			}
			if(cont==0){
				auxL[j] = i;
				position = j;
				// set left boundaries, when the index j is the first child.
				while(position && (position+1)%2==0){
					father = (position+1)>>1;
					if ((position+1)%2 > 1)
						father++;
					father--;
					auxL[father] = i;
					position = father;
				}
				cont=1;
			}else{
				cont++;
				if(cont == BLK){
					auxR[j] = i;
					position = j;
					// set right boundaries, when the index j is the last child. The last child always has rest == 1.
					while(position && (position+1)%2==1){
						father = (position+1)>>1;
						if ((position+1)%2 > 1)
							father++;
						father--;
						auxR[father] = i;
						position = father;
					}
					node++;
					j++;
					cont = 0;
				}
			}
		}
		if (cont)
			auxR[node] = i;
	}

	if (false && TRACE && leaves>0){
		cerr << endl << "min-max Intervals..." << endl;
		for (i=0; i<cantN; i++)
			cerr << "[" << auxL[i] << "," << auxR[i] << "] ";
		cerr << endl;
	}

	int MIN_BCK = 0;
	if (leaves){
		ulong segment;
		long int currSumBck;
		Aux_BkM = new int[cantIN];
		MIN_BCK = 1;

		for(i=cantIN; i>0; i--){
			child = i<<1;		// child is the second child of node i
			if (child > cantIN){
				if (child > firstLeaf)
					node = child - firstLeaf;
				else
					node = child - cantIN + leavesBottom;
				currSumBck = 0;
				miniBck = 0; // this 'min' correspond to the maximum level in the tree
				for (j=0, rb=N8W64-1; j<N8BLK; j++){
					segment = (P[((node+1)*BLK-1-BST*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BST*rb);

					if (currSumBck + T_MIN_BCK[segment] > miniBck)
						miniBck = currSumBck + T_MIN_BCK[segment];
					currSumBck += T_SUM_BLOCK[segment];
					if (rb == 0) rb=N8W64-1;
					else rb--;
				}
				for (j=0, rb=N8W64-1; j<N8BLK; j++){
					segment = (P[(node*BLK-1-BST*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BST*rb);

					if (currSumBck + T_MIN_BCK[segment] > miniBck)
						miniBck = currSumBck + T_MIN_BCK[segment];
					currSumBck += T_SUM_BLOCK[segment];
					if (rb == 0) rb=N8W64-1;
					else rb--;
				}
			}else{
				miniBck = Aux_BkM[child];
				currSumBck = sumAtPos(auxR[child]) - sumAtPos(auxL[child]-1);
				if (currSumBck + Aux_BkM[child-1] > miniBck)
					miniBck = currSumBck + Aux_BkM[child-1];

			}
			Aux_BkM[i-1] = miniBck;
			if (miniBck > MIN_BCK)
				MIN_BCK = miniBck;

		}
		//cerr << "Deleting auxL[] and auxR[]..." << endl;
		delete [] auxL;
		delete [] auxR;
	}

	lgBkM = 1 + (uint)(log(MIN_BCK)/log(2));
	if(cantIN){
		sizeDS = cantIN*lgBkM/W64;
		if ((cantIN*lgBkM)%W64)
			sizeDS++;
		//cerr << "Creating BkM[]... " << sizeDS*sizeof(ulong) << " Bytes" << endl;
		BkM = new ulong[sizeDS];
		sizeDS = sizeDS*sizeof(ulong);
		sizeRMM += sizeDS;

		ulong cMIN_BCK = 0; // count
		for(i=0; i<cantIN; i++){
			setNum64(BkM, cMIN_BCK, lgBkM, Aux_BkM[i]);
			cMIN_BCK += lgBkM;
		}
		//cerr << "Deleting Aux_BkM[]..." << endl;
		delete [] Aux_BkM;
	}else
		lgBkM=sizeDS=0;

	if (TRACE || SHOW_SIZE)
		cerr << " ** size of BkM[] " << sizeDS << " Bytes" << endl;
	cerr << " ** Total RMQRMM64 size: " << sizeRMM << " Bytes = " << (float)sizeRMM/(1024.0*1024.0) << " MB." << endl;
	

	if (TRACE){
		cerr << "MIN_BCK " << MIN_BCK << ", lgBkM " << lgBkM << endl;
		if (leaves>0) printTree();
		cerr << "# blocks " << nb << ", # blocks in the tree " << nb << endl;
		cerr << "nBin: " << nBin << ", Length last block " << lenLast << endl;
		cerr << "Sum relative to binary tree: " << sumAtPos(nBin) << endl;
		cerr << "======================"<< endl;
	}

	if (RUNTEST){
		test_search_min_block();
		test_sumAtPos();
		test_rank_1();
		test_select_1();
		test_positionMinblock();
		if (leaves>0) test_rmqi();
	}
}

void RMQRMM64::createTables(){
	ulong sizeDS, MAX_B, MAX_SumB, semiSum;
	ulong i, j, jSS, cont, rb, segment, posBLK, rOnes, nextSa;
	long int semiSumBlock, Min, auxMin;

	nBLK = nP/BLK;
	//cerr << "nBLK " << nBLK << endl;
	lenSS = (nP/2)/SS+1;
	//cerr << "lenSS " << lenSS << endl;
	posBLK = semiSum = MAX_B = MAX_SumB = 0;

	// here.... always sum >= 0, because is a BP sequence
	for (i=0; i<nBLK; i++){
		semiSumBlock = Min = 0;
		rb=N8W64-1;
		for (j=N8BLK; j>0; j--){
			segment = (P[(i*BLK+BST*(j-1))/W64] & RMMMasks[rb]) >> (W64m8-BST*rb);
			//printBitsUlong(segment);cerr<<endl;
			auxMin = semiSumBlock + T_MIN_BCK[segment];
			if(auxMin > Min)
				Min = auxMin;

			semiSumBlock += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		semiSumBlock >>= 1;
		semiSum += semiSumBlock;
		if (semiSum > MAX_SumB)
			MAX_SumB = semiSum;
		if (Min > (int)MAX_B)
			MAX_B = Min;
		posBLK += BLK;
	}

	MAX_B++;
	lg_MinB = 1 + (uint)(log(MAX_B)/log(2));
	lg_SumB = 1 + (uint)(log(MAX_SumB)/log(2));
	lg_SS = 1 + (uint)(log(nBLK)/log(2));

	if (TRACE)
		cerr << "MAX_B " << MAX_B << ", lg_MinB " << lg_MinB  << ", MAX_SumB " << MAX_SumB << ", lg_SumB " << lg_SumB << ", lg_SS " << lg_SS << endl;

	// W64 = 64 bits, and sizeof(ulong) = 8 bytes = 64 bits too
	cont = (nBLK+1)*lg_SumB/W64;
	if (((nBLK+1)*lg_SumB)%W64)
		cont++;
	TSumB = new ulong[cont];
	sizeDS = cont*sizeof(ulong);
	sizeRMM += sizeDS;
	if (TRACE || SHOW_SIZE) cerr << " ** size of TSumB[] " <<  sizeDS << " Bytes" << endl;
	setNum64(TSumB, 0, lg_SumB, 0);

	if(nBLK){
		cont = nBLK*lg_MinB/W64;
		if ((nBLK*lg_MinB)%W64)
			cont++;
		TMinB = new ulong[cont];
		sizeDS = cont*sizeof(ulong);
		sizeRMM += sizeDS;
		if (TRACE || SHOW_SIZE) cerr << " ** size of TMinB[] " <<  sizeDS << " Bytes" << endl;
	}else
		if (TRACE || SHOW_SIZE) cerr << " ** size of TMinB[] 0 Bytes" << endl;

	if(lenSS){
		cont = lenSS*lg_SS/W64;
		if ((lenSS*lg_SS)%W64)
			cont++;

		TSS = new ulong[cont];
		cont *= sizeof(ulong);
		sizeDS = cont;
		sizeRMM += sizeDS;
		if (TRACE || SHOW_SIZE) cerr << " ** size of TSS[] " << sizeDS << " Bytes" << endl;
	}else
		if (TRACE || SHOW_SIZE) cerr << " ** size of TSS[] 0 Bytes" << endl;

	posBLK = semiSum = 0;
	nextSa = SS;
	jSS = 1;
	for (i=0; i<nBLK; i++){
		semiSumBlock = Min = 0;
		rb=N8W64-1;
		for (j=N8BLK; j>0; j--){
			segment = (P[(i*BLK+BST*(j-1))/W64] & RMMMasks[rb]) >> (W64m8-BST*rb);
			auxMin = semiSumBlock + T_MIN_BCK[segment];
			if(auxMin > Min)
				Min = auxMin;
			semiSumBlock += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		semiSumBlock >>= 1;
		semiSum += semiSumBlock;
		setNum64(TSumB, (i+1)*lg_SumB, lg_SumB, semiSum);
		setNum64(TMinB, i*lg_MinB, lg_MinB, Min);

		posBLK += BLK;
		rOnes = semiSum + (posBLK/2);
		//cerr << "rOnes " << rOnes << " SB " << jS << " position " << posMin << endl;
		while (jSS*SS <= rOnes && jSS < (ulong)lenSS){
			setNum64(TSS, jSS*lg_SS, lg_SS, i);
			nextSa += SS;
			jSS++;
		}
	}
	while (jSS < (ulong)lenSS){
		setNum64(TSS, jSS*lg_SS, lg_SS, i);
		jSS++;
	}

	if (TRACE){
		cerr << endl << "TMinB[1.." <<nBLK<< "]... ";
		for (i=0; i<nBLK; i++)
			cerr << getNum64(TMinB, i*lg_MinB, lg_MinB) << " ";
		cerr << endl;
		cerr << "TSumB[1.." <<nBLK+1<< "]... ";
		for (i=0; i<=nBLK; i++)
			cerr << getNum64(TSumB, i*lg_SumB, lg_SumB) << " ";
		cerr << endl;
		cerr << "TSS[1.." <<lenSS<< "]... ";
		for (i=0; i<lenSS; i++)
			cerr << getNum64(TSS, i*lg_SS, lg_SS) << " ";
		cerr << endl << endl;
	}
}

// *******************************************************************************
// ********************************* BASIC OPERATIONS ****************************


// give the excess from 0 to pos
long int RMQRMM64::sumAtPos(long int pos){
	ulong rb, q;
	ulong blk = pos>>BBLK;
	long int j, l, sum=getNum64(TSumB, blk*lg_SumB, lg_SumB)<<1;

	j = blk<<BBLK;
	l=pos+1-BST;
	rb=0;
	while (j <= l){
		q = (P[j>>BW64] & RMMMasks[rb]) >> (W64m8-BST*rb);
		sum += T_SUM_BLOCK[q];
		if (rb == N8W64-1){
			rb=0;
			blk++;
		}
		else rb++;
		j +=BST;
	}
	while (j <= pos){
		if (readBit64(P, j))
			sum++;
		else
			sum--;
		j++;
	}

	return sum;
}

void RMQRMM64::test_sumAtPos(){
	ulong k;
	long int sum, sumPos;

	cerr << "RMQRMM64::test_sumAtPos..." << endl;
	/*k=128;
	for (ulong j=sum=0; j<=k; j++){
		if(readBit64(P, j))
			sum++;
		else
			sum--;
	}
	cerr << "sumAtPos(" << k << ") = " << sum << endl;
	sumPos = sumAtPos(k);
	if (sum != sumPos){
		cerr << "ERROR !! sumAtPos(" << k << ") = " << sumPos << " != sum = " << sum << endl;
		exit(1);
	}*/

	sum=0;
	for (k=0; k<(nP-2); k++){
		if(readBit64(P, k))
			sum++;
		else
			sum--;

		sumPos = sumAtPos(k);
		if (sum != sumPos){
			cerr << "ERROR !! sumAtPos(" << k << ") = " << sumPos << " != sum = " << sum << endl;
			exit(1);
		}
	}
	cerr << "  test_sumAtPos OK !!" << endl;
}



// for 0 <= i < n
ulong RMQRMM64::rank_1(ulong i){
	if(i >= nP-1) return nP>>1;
	ulong b, rest, q;
	ulong blk=(i+1)>>BBLK;
	ulong rank = (blk<<BMBLK) + getNum64(TSumB, blk*lg_SumB, lg_SumB);

	ulong x = blk<<BBLK;
	while(x+BSTMOne <= i){
		b = x>>BW64;
		rest = ((x+BST)%W64);
		q = (P[b] >> (W64-rest)) & 0xff;
		rank += __popcount_tab[q];
		x += BST;
	}

	// check last segment (< S) bit by bit...
	while (x<=i){
		if(readBit64(P,x))
			rank++;
		x++;
	}

	return rank;
}

void RMQRMM64::test_rank_1(){
	long int sumR, rank;
	ulong k;
	cerr << "RMQRMM64::test_rank_1..." << endl;

	/*k=1031;
	sumR=0;
	for (ulong j=0; j<=k; j++){
		if(readBit64(P, j))
			sumR++;
	}
	rank = rank_1(k);
	if (sumR != rank){
		cerr << "ERROR !! rank1(" << k << ") = " << rank << " != sumR = " << sumR << endl;
		exit(1);
	}*/

	sumR=0;
	for (k=0; k<nP; k++){
		if(readBit64(P, k))
			sumR++;

		rank = rank_1(k);
		if (sumR != rank){
			cerr << "ERROR !! rank1(" << k << ") = " << rank << " != sumR = " << sumR << endl;
			exit(1);
		}
	}
	cerr << "  test_rank_1 OK !!" << endl;
}


ulong RMQRMM64::select_1(ulong i){
	ulong j, jj, b, s, r, rr, sum, q;

	b = i>>BSS;
	s = getNum64(TSS, b*lg_SS, lg_SS);

	// Accumulate sum from super blocks
	r = sum = j = 0;
	if(s){
		b = s*lg_SumB;
		sum = getNum64(TSumB, b, lg_SumB);
		j = s<<BBLK;
		r = sum + (j>>1);

		// next SB
		b += lg_SumB;
		jj = j + BLK;
		rr = getNum64(TSumB, b, lg_SumB)+(jj>>1);
		while(rr < i){
			r = rr;
			j = jj;
			s++;

			b += lg_SumB;
			jj += BLK;
			rr = getNum64(TSumB, b, lg_SumB)+(jj>>1);
		}
	}

	// Accumulate sum block by block
	b = j>>BW64;
	s = (j+BST)%W64;
	q = (P[b] >> (W64-s)) & 0xff;
	rr = r+__popcount_tab[q];
	while(rr < i){
		r = rr;
		j += BST;

		b = j>>BW64;
		s = (j+BST)%W64;
		q = (P[b] >> (W64-s)) & 0xff;
		rr += __popcount_tab[q];
	}

	// check last segment (< S) bit by bit...
	while (r<i){
		if(readBit64(P,j))
			r++;
		j++;
	}

	return j-1;
}

void RMQRMM64::test_select_1(){
	ulong k, sumS, pos;

	cerr << "RMQRMM64::test_select_1..." << endl;
	/*k=32834;
	ulong j=0;
	for (sumS=0; sumS<k; j++){
		if(readBit64(P, j))
			sumS++;
	}
	j--;
	pos = select_1(k);
	cerr << "select_1(" << k <<") = " << pos << ", ? j=" << j << endl;
	if (pos != j){
		cerr << "ERROR !! select1(" << sumS << ") = " << pos << " != j = " << j << endl;
		exit(1);
	}
	exit(0);*/

	sumS=0;
	for (k=0; k<nP; k++){
		if(readBit64(P, k)){
			sumS++;

			pos = select_1(sumS);
			if (pos != k){
				cerr << "ERROR !! select1(" << sumS << ") = " << pos << " != k = " << k << endl;
				exit(1);
			}
		}
	}
	cerr << "  test_select_1 OK !!" << endl;
}


// give the excess from 0 to pos
long int RMQRMM64::sumAtBlock(long int blk){

	return (getNum64(TSumB, blk*lg_SumB, lg_SumB)<<1);

}

// return the position in the block 'blk' where is the minimum 'Min' of the block
ulong RMQRMM64::positionMinblock(ulong blk){
	int sum, min, Min;
	uint rest, q, rb;
	ulong b, posMin, x;

	Min = getNum64(TMinB, blk*lg_MinB, lg_MinB);
	posMin = x = ((blk+1)<<BBLK)-1;
	sum = rest = 0;
	rb = 1;

	b = x>>BW64;
	q = (P[b] >> rest) & 0xff;
	min = T_MIN_BCK[q];
	while ((min+sum) != Min){
		sum += T_SUM_BLOCK[q];
		rest += BST;
		x -= BST;
		if (rb == N8W64){
			rb=1;
			b--;
			rest = 0;
		}else
			rb++;

		q = (P[b] >> rest) & 0xff;
		min = T_MIN_BCK[q];
	}

	posMin = x-T_BCK_D[q][min-1];

	return posMin;
}

void RMQRMM64::test_positionMinblock(){
	long int sum, Min, k, j;
	ulong posMin;

	cerr << "RMQRMM64::test_positionMinblock()..." << endl;
	for (ulong t=0; t<leaves; t++){
		if (getNum64(TMinB, t*lg_MinB, lg_MinB) <= 0)
			continue;

		j = t*BLK;
		k = (t+1)*BLK-1;

		for (sum = 0;k>=j && !readBit64(P, k); k--, sum--);
		sum++;
		Min = sum;
		posMin = k;

		if (k>=j){
			for (k--; k>=j; k--){
				if(readBit64(P, k)){
					sum++;
					if (sum > Min){
						Min = sum;
						posMin = k;
					}
				}else
					sum--;
			}

			ulong position = positionMinblock(t);
			if (position != posMin){
				cerr << "ERROR !!... positionMinblock(" << t << ") = " << position << " != " << posMin << endl;
				cerr << "minimum = " << Min << endl;
				cerr << "TMinB[t] = " << getNum64(TMinB, t*lg_MinB, lg_MinB) << endl;
				exit(1);
			}
		}
	}

	cerr << "  positionMinblock() OK !!" << endl;
}


// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
long int RMQRMM64::computeSumOfNode(ulong node, ulong dist){
	long int sum = 0;
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	//sum += sumAtBlock(end) - sumAtBlock(ini);
	sum += (getNum64(TSumB, end*lg_SumB, lg_SumB)-getNum64(TSumB, ini*lg_SumB, lg_SumB))<<1;

	return sum;
}
// give the number of leaves of this node 'node=preorder+1' that has a distance 'dist' to the tree's depth
ulong RMQRMM64::leavesOfNode(ulong node, ulong dist, ulong *leafL){
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	*leafL = ini;
	return end-ini;
}

// give the number of leaves of this node 'node=preorder+1' that has a distance 'dist' to the tree's depth
ulong RMQRMM64::computeLeavesOfNode(ulong node, ulong dist){
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	return end-ini;
}

// search the rightmost minimum from x2 to x1 (sequential search)
void RMQRMM64::search_min_block(ulong x1, ulong x2, long int *min, long int *curSum, ulong *position){
	int Min, sum, auxMin;
	uint rest, q;
	ulong b, len = x2-x1+1, posMin = *position;

	Min = *min;
	sum = *curSum;
	while(len > BSTMOne){
		b = x2>>BW64;
		rest = (x2+1)%W64;
		if (b == (x2-BSTMOne)>>BW64)				// x and (x-7) are in the same word...
			q = (P[b] >> (W64-rest)) & 0xff;
		else
			q = ((P[b-1] << rest) | (P[b] >> (W64-rest))) & 0xff;
		auxMin = T_MIN_BCK[q];
		if (Min < sum+auxMin){
			Min = sum+auxMin;
			posMin = x2-T_BCK_D[q][auxMin-1];
		}
		sum += T_SUM_BLOCK[q];
		len -= BST;
		x2 -= BST;
	}
	if (len){
		// check last segment (len < 8) bit by bit...
		while (len>=1){
			if(readBit64(P,x2)){
				sum++;
				if (sum > Min){
					Min = sum;
					posMin = x2;
				}
			}else sum--;
			x2--;
			len--;
		}
	}
	*curSum = sum;
	*min = Min;
	*position = posMin;
}

void RMQRMM64::test_search_min_block(){
	long int sum, Min;
	ulong posMin, k, i, j;

	ulong TTEST = 1000;

	cerr << "RMQRMM64::test_search_min_block()..." << endl;
	for (ulong t=0; t<TTEST; t++){
		i = 1+(rand() % (nP/2));
		j = i+(rand()%(nP-i))-1;
		if (i>j) continue;

		for (k=j, sum=0; k>=i && !readBit64(P,k); k--, sum--);
		sum++;
		Min=sum;
		posMin=k;
		if (k<i) continue;

		for (k--; k>=i; k--){
			if(readBit64(P, k)){
				sum++;
				if (sum > Min){
					Min = sum;
					posMin = k;
				}
			}else
				sum--;
		}

		//if (TRACE) cerr << "search_min_block(" << i << ", " << j << ") " << endl;
		long int min, curSum;
		ulong position;
		min = curSum = 0;

		for (k=j, curSum=0; k>=i && !readBit64(P,k); k--, curSum--);
		curSum++;
		min=curSum;
		position=k;
		search_min_block(i,k-1,&min, &curSum, &position);
		if (position != posMin || min != Min){
			cerr << "ERROR !!... search_min_block(" << i << ", " << j << ") = " << position << " != " << posMin << endl;
			cerr << "minimum found = " << min << " =? " << Min << endl;
			exit(1);
		}
	}

	cerr << "  test_search_min_block() OK !!" << endl;
}

// return the position of the open parenthesis closet to the root between i and j
ulong RMQRMM64::rmqi_rmm(ulong x1, ulong x2, long int *min, long int *currSum, ulong posMin){
	ulong leafL, lLeaf, aLeaf, node, nodeMin, posM, g, dNode, dist = 1;
	long int minT, sumMin;
	bool isNode = false;
	bool isLeaf = false;

	lLeaf = x1>>BBLK;
	node = x2>>BBLK;

	// [1]- search the minimum inside the rightmost block...
	if (lLeaf == node){
		// here is the answer
		search_min_block(x1, x2, min, currSum, &posMin);
		return posMin;
	}else
		search_min_block((node<<BBLK), x2, min, currSum, &posMin);
	long int Min = *min, sum = *currSum;

	// [2]- We looking for in the next leaf
	if (node%2){	// this is a right leaf
		node--;
		minT = getNum64(TMinB, node*lg_MinB, lg_MinB)+sum;
		if (node == lLeaf){		// here is the answer
			if(minT > Min)
				search_min_block(x1, ((node+1)<<BBLK)-1, min, currSum, &posMin);
			return posMin;
		}

		if (minT > Min){
			nodeMin = node;
			sumMin = sum;
			isNode = isLeaf = true;
			Min = minT;
		}
		sum += (getNum64(TSumB, (node+1)*lg_SumB, lg_SumB)-getNum64(TSumB, node*lg_SumB, lg_SumB))<<1;
	}
	if (node == 0)
		return posMin;

	// [3]- We climb recomputing the Min until the position i.
	aLeaf = node;
	if (node < leavesBottom){
		if (node == leavesBottom-1)
			node = cantIN-1;
		else
			node += firstLeaf;
	}else
		node += cantIN-leavesBottom;
	node>>= 1;
	node--;
	dist++;
	g = leavesOfNode(node+1, dist, &leafL);
	while (lLeaf < aLeaf-g){                  /////// group > aLeaf && ...
		aLeaf -= g;
		minT = sum+getNum64(BkM, node*lgBkM, lgBkM);
		if (minT > Min){
			Min = minT;
			isNode = true;
			isLeaf = false;
			nodeMin = node;
			sumMin = sum;
			dNode = dist;
		}
		sum += (getNum64(TSumB, (leafL+g)*lg_SumB, lg_SumB)-getNum64(TSumB, aLeaf*lg_SumB, lg_SumB))<<1;
		if (node%2==0)
			node--;
		else{
			node>>=1;		// go to my uncle
			node--;
			dist++;
			g<<=1;
		}
		g = leavesOfNode(node+1, dist, &leafL);
	}

	// [4]- We move down recomputing the Min and reach to the left boundaries leaf.
	node++;
	node<<= 1;
	dist--;
	while (node < cantIN){
		g = leavesOfNode(node+1, dist, &leafL);
		if(lLeaf < aLeaf-g){
			aLeaf -= g;
			minT = sum+getNum64(BkM, node*lgBkM, lgBkM);
			if (minT > Min){
				Min = minT;
				isNode = true;
				isLeaf = false;
				nodeMin = node;
				sumMin = sum;
				dNode = dist;
			}
			sum += (getNum64(TSumB, (leafL+g)*lg_SumB, lg_SumB)-getNum64(TSumB, aLeaf*lg_SumB, lg_SumB))<<1;
			node--;
		}else{
			node++;
			node<<= 1;
			dist--;
		}
	}
	if (node > firstLeaf)			// at this point, 'node' is a leaf; node >= firstLeaf
		node -= firstLeaf;
	else
		node += leavesBottom - cantIN;

	// [5] looking for in the leaves of the leftmost side
	if (node == lLeaf){
		minT = getNum64(TMinB, node*lg_MinB, lg_MinB)+sum; // ESTA BIEN node ??
		if (isNode){
			if (minT > Min){
				// puede ser que el minimo este a la izquierda de x1 en esta hoja, validar !!
				posM=posMin;
				search_min_block(x1, ((node+1)<<BBLK)-1, &Min, &sum, &posM);
				if (posM < posMin)
					return posM;
			}
		}else{
			if (minT > Min)
				search_min_block(x1, ((node+1)<<BBLK)-1, &Min, &sum, &posMin);
			return posMin;
		}
	}else{
		minT = getNum64(TMinB, node*lg_MinB, lg_MinB)+sum;

		if (isNode){
			if (minT > Min){
				Min = minT;
				posMin = positionMinblock(node);
				sum += (getNum64(TSumB, (node+1)*lg_SumB, lg_SumB)-getNum64(TSumB, node*lg_SumB, lg_SumB))<<1;
				search_min_block(x1, (node<<BBLK)-1, &Min, &sum, &posMin);
				return posMin;
			}
			// node is a right leaf
			sum += (getNum64(TSumB, (node+1)*lg_SumB, lg_SumB)-getNum64(TSumB, node*lg_SumB, lg_SumB))<<1;
			minT = getNum64(TMinB, (node-1)*lg_MinB, lg_MinB)+sum;

			if (minT > Min){
				posM=posMin;
				search_min_block(x1, (node<<BBLK)-1, &Min, &sum, &posM);
				if (posM < posMin)
					return posM;
			}
		}else{
			if (minT > Min){
				isLeaf = true;
				nodeMin = node;
				sumMin = sum;
				Min = minT;
			}

			sum += (getNum64(TSumB, (node+1)*lg_SumB, lg_SumB)-getNum64(TSumB, node*lg_SumB, lg_SumB))<<1;
			minT = getNum64(TMinB, (node-1)*lg_MinB, lg_MinB)+sum;
			if (minT > Min){
				// puede ser que el minimo este a la izquierda de x1 en esta hoja, validar !!
				posM=posMin;
				search_min_block(x1, (node<<BBLK)-1, &Min, &sum, &posM);
				if (posM < posMin)
					return posM;
			}
			if (isLeaf)
				posMin = positionMinblock(nodeMin);

			return posMin;
		}
	}

	if(isLeaf)
		return positionMinblock(nodeMin);

	// [6] Here the minimum value is in a block which descent of the internal node 'nodeMin',
	// then we must descent in order to find its position.
	nodeMin = (nodeMin+1)<<1;
	dNode--;
	sum = sumMin;
	while (nodeMin < cantIN){
		minT = getNum64(BkM, nodeMin*lgBkM, lgBkM);
		sumMin = computeSumOfNode(nodeMin+1, dNode);
		if (Min == sum+minT){   // is it ? then go to the right without updating the Min
			nodeMin++;
			nodeMin<<= 1;
			dNode--;
		}else{					// else, go to the left updating the sum
			nodeMin--;
			sum += sumMin;
		}
	}
	if (nodeMin > firstLeaf)
		nodeMin -= firstLeaf;
	else
		nodeMin += leavesBottom - cantIN;

	// Here, I get the right leaf... nodeMin is a right leaf
	sumMin = (getNum64(TSumB, (nodeMin+1)*lg_SumB, lg_SumB)-getNum64(TSumB, nodeMin*lg_SumB, lg_SumB))<<1;
	sumMin += getNum64(TMinB, (nodeMin-1)*lg_MinB, lg_MinB);

	if (sumMin > (long int)getNum64(TMinB, nodeMin*lg_MinB, lg_MinB))
		return positionMinblock(nodeMin-1);
	return positionMinblock(nodeMin);

}

// always return a position of open parenthesis with the minimum value
ulong RMQRMM64::rmqi(ulong i, ulong j){
	long int min, sum;
	min = sum = 0;
	if(j < nBin)
		return rmqi_rmm(i, j, &min, &sum, j);
	else{
		ulong posMin = j;

		if(i >= nBin){
			search_min_block(i, j, &min, &sum, &posMin);
			return posMin;
		}else{
			if (nBin){
				search_min_block(nBin, j, &min, &sum, &posMin);
				return rmqi_rmm(i, nBin-1, &min, &sum, posMin);
			}
		}
	}
	return 0;
}

void RMQRMM64::test_rmqi(){
	ulong i, j, posMin, posRmqi;
	long int sum, Min, k;

	/*i=65400, j=65598;
	for (sum=0, Min=0, posMin=i, k=j; k>=(int)i; k--){
		if(readBit64(P, k)){
			sum++;
			if (sum>Min){
				Min = sum;
				posMin = k;
			}
		}else
			sum--;
	}
	cerr << "MINIMUM(" << i << ", " << j << ") = " << posMin << endl;
	posRmqi = rmqi(i,j);
	if (posRmqi != posMin){
		cerr << "ERROR !!... rmqi(" << i << ", " << j << ") = " << posRmqi << " != " << posMin << endl;
		exit(1);
	}
	exit(0);*/


	cerr << "RMQRMM64::test_rmqi..." << endl;
	i=0;
	ulong sumTest = 1;
	if(nP >= 10000)
		sumTest = 1000;
	if(nP >= 1000000)
		sumTest = 10000;
	if(nP >= 10000000)
		sumTest = 1000000;
	if(nP >= 100000000)
		sumTest = 10000000;
	while(i<nP){
		while(i<nP && readBit64(P, i)==0)
			i++;
		j=i+1;
		while(j<nP){
			while(j<nP && readBit64(P, j)==0)
				j++;

			if (j<nP){
				for (sum=0, Min=0, posMin=i, k=j; k>=(long int)i; k--){
					if(readBit64(P, k)){
						sum++;
						if (Min < sum){
							Min = sum;
							posMin = k;
						}
					}else
						sum--;
				}

				//if (TRACE) cerr << "rmqi(" << i << ", " << j << ") " << endl;
				posRmqi = rmqi(i,j);
				if (posRmqi != posMin){
					cerr << "ERROR !!... rmqi(" << i << ", " << j << ") = " << posRmqi << " != " << posMin << endl;
					exit(1);
				}
				//if (TRACE) cerr << "= " << posRmqi << endl;
			}
			j+=sumTest;
		}
		i+=sumTest;
	}

	cerr << "test_rmqi OK !!" << endl;
}

// =================================================================================================================

ulong RMQRMM64::queryRMQ(ulong i, ulong j){
	return rank_1(rmqi(select_1(i+2),select_1(j+2)))-2;
}

uint RMQRMM64::getSize(){
	return sizeRMM;
}

void RMQRMM64::saveDS(char *fileName){
	cerr << "Save data structure in " << fileName << endl;
	ofstream os (fileName, ios::binary);
	cerr << "   Data structure size: " << sizeRMM << endl;

	if(TRACE){
		cerr << "Variables load: " << endl;
		cerr << "nP " << nP << endl;
		cerr << "nW " << nW << endl;
		cerr << "nBin " << nBin << endl;
		cerr << "cantN " << cantN << endl;
		cerr << "cantIN " << cantIN << endl;
		cerr << "leaves " << leaves << endl;
		cerr << "leavesBottom " << leavesBottom << endl;
		cerr << "firstLeaf " << firstLeaf << endl;

		cerr << "lenLB " << lenLast << endl;
		cerr << "nBLK " << nBLK << endl;
		cerr << "lenSS " << lenSS << endl;
		cerr << "lg_MinB " << lg_MinB << endl;
		cerr << "lg_SumB " << lg_SumB << endl;
		cerr << "lg_SS " << lg_SS << endl;
		cerr << "lgBkM " << lgBkM << endl;
	}

	os.write((const char*)&nP, sizeof(ulong));
	os.write((const char*)&nW, sizeof(ulong));
	os.write((const char*)&nBin, sizeof(ulong));
	os.write((const char*)&cantN, sizeof(ulong));
	os.write((const char*)&cantIN, sizeof(ulong));
	os.write((const char*)&leaves, sizeof(ulong));
	os.write((const char*)&leavesBottom, sizeof(ulong));
	os.write((const char*)&firstLeaf, sizeof(ulong));
	os.write((const char*)&nBLK, sizeof(ulong));

	os.write((const char*)&lenLast, sizeof(uint));
	os.write((const char*)&lenSS, sizeof(uint));
	os.write((const char*)&lg_MinB, sizeof(uint));
	os.write((const char*)&lg_SumB, sizeof(uint));
	os.write((const char*)&lg_SS, sizeof(uint));
	os.write((const char*)&lgBkM, sizeof(uint));

	ulong sizeDT = 3072; // =  2*512 + 2048;  size for T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[]
	if(TRACE) cerr << " .- T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[] " << sizeDT << " Bytes" << endl;

	ulong size = nP >> BW64;
	if (nP % W64)
		size++;
	os.write((const char*)P, size*sizeof(ulong));				// save P[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cerr << " .- P[] " << size*sizeof(ulong) << " Bytes" << endl;

	if(cantIN){
		size = cantIN*lgBkM/W64;
		if ((cantIN*lgBkM)%W64)
			size++;
		os.write((const char*)BkM, size*sizeof(ulong));			// save BkM[]
		sizeDT += size*sizeof(ulong);
		if(TRACE) cerr << " .- BkM[] " << size*sizeof(ulong) << " Bytes" << endl;
	}else
		if(TRACE) cerr << " .- BkM[] 0 Bytes" << endl;

	size = (nBLK+1)*lg_SumB/W64;
	if (((nBLK+1)*lg_SumB)%W64)
		size++;
	os.write((const char*)TSumB, size*sizeof(ulong));			// save TSumB[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cerr << " .- TSumB[] " << size*sizeof(ulong) << " Bytes" << endl;

	if(lenSS){
		size = lenSS*lg_SS/W64;
		if ((lenSS*lg_SS)%W64)
			size++;
		os.write((const char*)TSS, size*sizeof(ulong));			// save TSS[]
		sizeDT += size*sizeof(ulong);
		if(TRACE) cerr << " .- TSS[] " << size*sizeof(ulong) << " Bytes" << endl;
	}else
		if(TRACE) cerr << " .- TSS[] 0 Bytes" << endl;

	if(nBLK && leaves){
		size = nBLK*lg_MinB/W64;
		if ((nBLK*lg_MinB)%W64)
			size++;
		os.write((const char*)TMinB, size*sizeof(ulong));		// save TMinB[]
		sizeDT += size*sizeof(ulong);
		if(TRACE) cerr << " .- TMinB[] " << size*sizeof(ulong) << " Bytes" << endl;
	}else
		if(TRACE) cerr << " .- TMinB[] 0 Bytes" << endl;

	os.close();
	cerr << "   Total bytes saved from data structure: " << sizeDT << endl;
}

void RMQRMM64::loadDS(char *fileName){
	cerr << " Load data structure from " << fileName << endl;
	ifstream is(fileName, ios::binary);
	cout << " Load data structure from " << fileName << endl;

	is.read((char*)&nP, sizeof(ulong));
	is.read((char*)&nW, sizeof(ulong));
	is.read((char*)&nBin, sizeof(ulong));
	is.read((char*)&cantN, sizeof(ulong));
	is.read((char*)&cantIN, sizeof(ulong));
	is.read((char*)&leaves, sizeof(ulong));
	is.read((char*)&leavesBottom, sizeof(ulong));
	is.read((char*)&firstLeaf, sizeof(ulong));
	is.read((char*)&nBLK, sizeof(ulong));

	is.read((char*)&lenLast, sizeof(uint));
	is.read((char*)&lenSS, sizeof(uint));
	is.read((char*)&lg_MinB, sizeof(uint));
	is.read((char*)&lg_SumB, sizeof(uint));
	is.read((char*)&lg_SS, sizeof(uint));
	is.read((char*)&lgBkM, sizeof(uint));

	if(TRACE){
		cerr << "Variables load: " << endl;
		cerr << "nP " << nP << endl;
		cerr << "nW " << nW << endl;
		cerr << "nBin " << nBin << endl;
		cerr << "cantN " << cantN << endl;
		cerr << "cantIN " << cantIN << endl;
		cerr << "leavesBottom " << leavesBottom << endl;
		cerr << "leaves " << leaves << endl;
		cerr << "firstLeaf " << firstLeaf << endl;

		cerr << "lenLB " << lenLast << endl;
		cerr << "nBLK " << nBLK << endl;
		cerr << "lenSS " << lenSS << endl;
		cerr << "lg_MinB " << lg_MinB << endl;
		cerr << "lg_SumB " << lg_SumB << endl;
		cerr << "lg_SS " << lg_SS << endl;
		cerr << "lgBkM " << lgBkM << endl;
	}

	sizeRMM = 3072; // = 2*512 + 2048;  size for T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[]
	if(TRACE) cerr << " .- size of T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[] = " << sizeRMM << " Bytes" << endl;

	ulong sizeDS = nP >> BW64;
	if (nP % W64)
		sizeDS++;
	P = new ulong[sizeDS];
	is.read((char*)P, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cerr << " .- P[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	if(cantIN){
		sizeDS = cantIN*lgBkM/W64;
		if ((cantIN*lgBkM)%W64)
			sizeDS++;
		BkM = new ulong[sizeDS];
		is.read((char*)BkM, sizeDS*sizeof(ulong));
		sizeRMM += sizeDS*sizeof(ulong);
		if(TRACE) cerr << " .- BkM[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;
	}else
		if(TRACE) cerr << " .- BkM[] 0 Bytes" << endl;

	sizeDS = (nBLK+1)*lg_SumB/W64;
	if (((nBLK+1)*lg_SumB)%W64)
		sizeDS++;
	TSumB = new ulong[sizeDS];
	is.read((char*)TSumB, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cerr << " .- TSumB[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	if(lenSS){
		sizeDS = lenSS*lg_SS/W64;
		if ((lenSS*lg_SS)%W64)
			sizeDS++;
		TSS = new ulong[sizeDS];
		is.read((char*)TSS, sizeDS*sizeof(ulong));
		sizeRMM += sizeDS*sizeof(ulong);
		if(TRACE) cerr << " .- TSS[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;
	}else
		if(TRACE) cerr << " .- TSS[] 0 Bytes" << endl;

	if(nBLK && leaves){
		sizeDS = nBLK*lg_MinB/W64;
		if ((nBLK*lg_MinB)%W64)
			sizeDS++;
		TMinB = new ulong[sizeDS];
		is.read((char*)TMinB, sizeDS*sizeof(ulong));
		sizeRMM += sizeDS*sizeof(ulong);
		if(TRACE) cerr << " .- TMinB[] " << sizeDS*sizeof(sizeDS) << " Bytes" << endl;
	}else
		if(TRACE) cerr << " .- TMinB[] 0 Bytes" << endl;

	is.close();
	cerr << " Data Structure loaded !!" << endl;
}

// *************************** PRINT THE TREE ************************************
void RMQRMM64::printTree(){
	ulong i, j, cant, acum, rb, segment;
	long int currSum, mini;
	uint cMIN_BCK;
	cMIN_BCK = 0;

	cerr << "Min-max Binary Tree...[min bck]i(sum)" << endl;
	cant = 1;
	for (acum=j=i=0; i<cantIN; i++, j++, acum++){
		if (j==cant){
			cerr << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cerr << " | ";
			acum = 0;
		}

		mini = getNum64(BkM, cMIN_BCK, lgBkM);
		cerr << i << "_[" << mini << "]";
		cMIN_BCK += lgBkM;
	}

	for (i=leavesBottom; i<leaves; i++, j++, acum++){
		if (j==cant){
			cerr << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cerr << " | ";
			acum = 0;
		}

		rb = 3;
		currSum = 0;
		mini = 0;
		for (j=0; j<N8BLK; j++){
			segment = (P[((i+1)*BLK-1-BST*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BST*rb);
			if (currSum + T_MIN_BCK[segment] > mini)
				mini = currSum + T_MIN_BCK[segment];
			currSum += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		cerr << i << "_h[" << mini << "](" << currSum << ") ";
	}
	cerr << endl;
	for (i=0; i<leavesBottom; i++, j++, acum++){
		if (j==cant){
			cerr << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cerr << " | ";
			acum = 0;
		}
		rb = 3;
		currSum = 0;
		mini = 0;
		for (j=0; j<N8BLK; j++){
			segment = (P[((i+1)*BLK-1-BST*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BST*rb);
			if (currSum + T_MIN_BCK[segment] > mini)
				mini = currSum + T_MIN_BCK[segment];
			currSum += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		cerr << i << "_h[" << mini << "]"<< i << "(" << currSum << ") ";
	}
	cerr << endl;
}

RMQRMM64::~RMQRMM64() {
	delete [] TSumB;
	if (nP) delete [] P;
	if(cantIN) delete [] BkM;
	if(lenSS) delete [] TSS;
	if(nBLK) delete [] TMinB;

	cerr << " ~ RMQRMM64 destroyed !!" << endl;
}

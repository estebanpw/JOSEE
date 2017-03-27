#include "alignment_functions.h"


static int64_t PAM[5][5]={4,-3,-3,-3,-3,-3,4,-3,-3,-3,-3,-3,4,-3,-3,-3,-3,-3,4,-3,-3,-3,-3,-3,-3};
int valOfNucl(char c){
	if(c=='A') return 0;
	if(c=='C') return 1;
	if(c=='G') return 2;
	if(c=='T') return 3;
	return 4;
}
/*
Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/

struct cell NWscore2rows(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell * mc, struct cell * f0, struct cell * f1){
    
    uint64_t i, j, k, iCounter=0, jCounter, currEgapR, currEgapG;
	int64_t scoreDiagonal = INT64_MIN, scoreLeft = INT64_MIN, scoreRight = INT64_MIN, score = INT64_MIN;
	struct cell * faux;
	
    struct cell mf;
    

    if(mc == NULL || f0 == NULL || f1 == NULL){
    	printf("Could not allocate memory.\n");
    	exit(-1);
    }
    
    
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");

    f0[iCounter].score = PAM[valOfNucl(X[Xstart])][valOfNucl(Y[Ystart])];
    f0[iCounter].igaps = 0;
    f0[iCounter].egaps = 0;
    f0[iCounter].ident = (((f0[iCounter].score) > (0)) ? (1) : (0));
    f0[iCounter].xs = Xstart;
    f0[iCounter].ys = Ystart;
    f0[iCounter].xe = Xstart;
    f0[iCounter].ye = Ystart;
    mc[iCounter] = f0[iCounter];
    //printf("    %03"PRId64" ", f0[iCounter].score);
    iCounter++;

    for(i=Ystart+1;i<Yend;i++){

        f0[iCounter].score = PAM[valOfNucl(X[Xstart])][valOfNucl(Y[i])] + iGap + ((i-Ystart)-1)*eGap;
    	f0[iCounter].igaps = 1;
    	f0[iCounter].egaps = i - Ystart;
    	f0[iCounter].ident = (((f0[iCounter].score) > (0)) ? (1) : (0));
    	f0[iCounter].xs = Xstart;
    	f0[iCounter].ys = i;
    	f0[iCounter].xe = Xstart;
    	f0[iCounter].ye = i;

    	//Set every column max
    	mc[iCounter] = f0[iCounter];
    	//goodPrint(f0[iCounter].score);
        //printf("%03"PRId64" ", f0[iCounter].score);
    	
    	iCounter++;
	}
	
	//Set row max
	mf = f0[0];
    //printf("\n");
	
	iCounter=1;

	//Go through full matrix with 2 rows
	for(i=Xstart+1;i<Xend;i++){
		//Fill first rowcell
		//printf("ROW: %"PRId64" ||",i);
		//printf("%.2f \n", ((float)i/Xend));
		
        f1[0].score = PAM[valOfNucl(X[i])][valOfNucl(Y[Ystart])] + iGap + ((i-Xstart)-1)*eGap;
		f1[0].xs = i;
		f1[0].ys = Ystart;
		f1[0].ident = (((PAM[valOfNucl(X[i])][valOfNucl(Y[Ystart])]) > (0)) ? (1) : (0));
        f1[0].igaps = 1;
        f1[0].egaps = i - (Xstart + 1);		
		f1[0].xe = i;
		f1[0].ye = Ystart;

		mf = f0[0];

		//goodPrint(f1[0].score);
        //printf("%03"PRId64" ", f1[0].score);
		jCounter=1;
		for(j=Ystart+1;j<Yend;j++){
		//for(j=redir[i][0]+1;j<redir[i][1];j++){
			//Check if max in row has changed
			if(jCounter > 1 && mf.score <= f0[jCounter-2].score){
				mf = f0[jCounter-2];
				mf.xe = i-1;
				mf.ye = j-2;
			}
			
			score = PAM[valOfNucl(X[i])][valOfNucl(Y[j])];
			scoreDiagonal = f0[jCounter-1].score + score;
			if(jCounter>1 && iCounter>=1){
				scoreLeft = mf.score + iGap + (j - (mf.ye+2))*eGap + score;
				currEgapR = (j - (mf.ye+2));
				}else{
					scoreLeft = INT64_MIN;
				}
				
			if(iCounter>1 && jCounter>=1){
				scoreRight = mc[jCounter-1].score + iGap + (i - (mc[j-1].xe+2))*eGap + score;
				currEgapG =  (i - (mc[j-1].xe+2));
				}else{
					scoreRight = INT64_MIN;
				}
			
			//Choose maximum
			//f1[jCounter] = max(max(scoreDiagonal,scoreLeft),scoreRight);
			
			if(scoreDiagonal >= MAX(scoreLeft, scoreRight)){
				//Diagonal
				f1[jCounter] = f0[jCounter-1];
				f1[jCounter].score = scoreDiagonal;
				if(PAM[valOfNucl(X[i])][valOfNucl(Y[j])] > 0) f1[jCounter].ident += 1;
								
			}else if(scoreRight >= scoreLeft){
				//Gap in genome
				f1[jCounter] = mc[jCounter-1];
				f1[jCounter].score = scoreRight;
				f1[jCounter].igaps += 1;
				f1[jCounter].egaps += currEgapG;
				
			}else{
				//Gap in read
				f1[jCounter] = mf;
				f1[jCounter].score = scoreLeft;
				f1[jCounter].igaps += 1;
				f1[jCounter].egaps += currEgapR;
			}

            if(iCounter>1 && jCounter>=1 && f0[jCounter].score >= mc[jCounter-1].score){
                mc[jCounter-1].score = f0[jCounter].score;
                mc[jCounter-1].xe = i-2;
                mc[jCounter-1].ye = j-1;
            }




			//Update movement
			f1[jCounter].xe = i;
			f1[jCounter].ye = j;
			//goodPrint(f1[jCounter].score);
            //printf("%03"PRId64" ", f1[jCounter].score);
			jCounter++;
		}
        //printf("\n");
        /*
		kCounter=0;
		for(k=Ystart;k<Yend;k++){
			//Update column maximum at j
			if(mc[kCounter].score <= f0[kCounter].score){
				mc[kCounter] = f0[kCounter];
				mc[kCounter].xe = i-1;
				mc[kCounter].ye = kCounter;
			}

			kCounter++;
		}
        */
		//Switch rows
		
		faux = f0;
		f0 = f1;
		f1 = faux;
		
		iCounter++;
	}



    int64_t bestScore=INT64_MIN, bestId = 0;
    for(k=Ystart;k<Yend;k++){
    	//printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", f0[k].xs, f0[k].ys, f0[k].xe, f0[k].ye, f0[k].score, f0[k].gaps, f0[k].ident);
    	if(f0[k].score >= bestScore){
    		bestScore = f0[k].score;
    		//printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", f0[k].xs, f0[k].ys, f0[k].xe, f0[k].ye, f0[k].score, f0[k].gaps, f0[k].ident);
    		bestId = k;
    	}
    }
    return f0[bestId];
}
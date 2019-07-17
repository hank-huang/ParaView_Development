//-----------------------------------------------------------------------------------//
//  Author and Copyright Owner of this Code: Mirza Faisal Beg                        //
//  Ph.D. Candidate in Biomedical Engineering, The Johns Hopkins University          //
//  School of Medicine. mfbeg@bme.jhu.edu                                            //
//  Code written towards thesis: Computing Metrics on Medical Images via Metrics     //
//  on Flows of Diffeomorphisms using Image Intensity and Landmark data.             //
//-----------------------------------------------------------------------------------//


#ifndef _MF_BYTE_SWAP_H_
#define _MF_BYTE_SWAP_H_

class classByteSwap{

public:


	static void BYTE_SWAP_INT(unsigned char* pntr){
	

		//before swap:  b0 b1 b2 b3
		//after  swap:  b3 b2 b1 b0
 
		unsigned char b0, b1, b2, b3;

		b0 = *(pntr  );
		b1 = *(pntr+1);
		b2 = *(pntr+2);
		b3 = *(pntr+3);

		*pntr = b3;
		*(pntr+1) = b2;
		*(pntr+2) = b1;
		*(pntr+3) = b0;

	}//BYTE_SWAP_INT(unsigned char* pntr)
	

        static void BYTE_SWAP_INT(int* pntr){

		classByteSwap::BYTE_SWAP_INT( (unsigned char*) pntr);

        }//BYTE_SWAP_INT(int* pntr)









        static void BYTE_SWAP_FLOAT(unsigned char* pntr){
        

                //before swap:  b0 b1 b2 b3
                //after  swap:  b3 b2 b1 b0

                unsigned char b0, b1, b2, b3;

                b0 = *(pntr  );
                b1 = *(pntr+1);
                b2 = *(pntr+2);
                b3 = *(pntr+3);

                *pntr     = b3;
                *(pntr+1) = b2;
                *(pntr+2) = b1;
                *(pntr+3) = b0;

        }//BYTE_SWAP_FLOAT(unsigned char* pntr)


        static void BYTE_SWAP_FLOAT(float* pntr){

		classByteSwap::BYTE_SWAP_FLOAT( (unsigned char*) pntr);

	}////BYTE_SWAP_FLOAT(float* pntr)








	static void BYTE_SWAP_SHORT(unsigned char* pntr){

		//before swap:   b0 b1
		//after swap: 	 b1 b0

		unsigned char b0, b1;
		
		b0 = *(pntr  );
                b1 = *(pntr+1);

		*pntr     = b1;
                *(pntr+1) = b0;

	}//BYTE_SWAP_SHORT(unsigned char* pntr)




        static void BYTE_SWAP_SHORT(short* pntr){

                classByteSwap::BYTE_SWAP_SHORT( (unsigned char*) pntr);

        }////BYTE_SWAP_FLOAT(float* pntr)





	 static void BYTE_SWAP_DOUBLE(unsigned char* pntr){

		//before swap: b0 b1 b2 b3 b4 b5 b6 b7
		//after swap:  b7 b6 b5 b4 b3 b2 b1 b0 

		unsigned char b0, b1, b2, b3, b4, b5, b6, b7;

                b0 = *(pntr  );
                b1 = *(pntr+1);
                b2 = *(pntr+2);
                b3 = *(pntr+3);
                b4 = *(pntr+4);
                b5 = *(pntr+5);
                b6 = *(pntr+6);
                b7 = *(pntr+7);

                *pntr     = b7;
                *(pntr+1) = b6;
                *(pntr+2) = b5;
                *(pntr+3) = b4;
                *(pntr+4) = b3;
                *(pntr+5) = b2;
                *(pntr+6) = b1;
                *(pntr+7) = b0;


	}//BYTE_SWAP_DOUBLE(unsigned char* pntr)



         static void BYTE_SWAP_DOUBLE(double* pntr){

			classByteSwap::BYTE_SWAP_DOUBLE((unsigned char*) pntr);

        }//BYTE_SWAP_DOUBLE(double* pntr)



	




};


#endif

/*------------------------------------------------------*/
/* Prog    : TpIFT6150-2.c                              */
/* Auteur  : Ahmed Amine DAOU & Wrushabh Warshe         */
/* Date    :                                            */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo2.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/                              
/*------------------------------------------------*/
#define NAME_IMG_IN  "photograph"
#define NAME_IMG_GAUSSIEN "TpIFT6150-2-gaussien"
#define NAME_IMG_GRADIENT "TpIFT6150-2-gradient"
#define NAME_IMG_SUPPRESSION "TpIFT6150-2-suppression"
#define NAME_IMG_CANNY "TpIFT6150-2-canny"
#define NAME_IMG_CANNY_AUTO "TpIFT6150-2-cannySemiAuto"
#define MAXIMUM(x, y) (((x) > (y)) ? (x) : (y))

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/                     
/*------------------------------------------------*/

/* 1-Rehaussement de Canny*/

  /* filtre gaussient
      params:
      ConvolutionR : matrice réelle de l'image flouée
      ConvolutionI : matrice imaginaire de l'image flouée

      imgINR  : matrice réelle a flouer
      imgINI  : matrice Imaginaire a flouer

      IMRGaussian: matrice réelle du filtre gaussien 
      IMRIaussian: matrice imaginaire du filtre gaussien

      sigma: ecart-type de la gaussienne 
      width & length: mesures des matrices
  */
void FiltreGaussien(float** ConvolutionR,float** ConvolutionI, float** imgINR,float** imgINI,float** IMRGaussian,float** IMIGaussian,float sigma,int length,int width){

  for (int i = 0; i < length; i++)
    {
      for (int j = 0; j < width; j++)
        {
          IMRGaussian[i][j]=funcgauss2D(i,j,sigma);
          IMIGaussian[i][j]=0.0;
          imgINI[i][j]=0.0;
        }
    }
    /* dans le domaine frequenciel*/
    FFTDD(imgINR, imgINI, length, width);
    FFTDD(IMRGaussian, IMRGaussian, length, width);

    /* application du filtre gaussien dans le domaine frequenciel*/ 
    MultMatrix(ConvolutionR, ConvolutionI, imgINR, imgINI, IMRGaussian, IMIGaussian, length, width);

    /*revenir dans le domaine spacial*/
    IFFTDD(ConvolutionR, ConvolutionI, length, width);
    SaveImagePgm(NAME_IMG_GAUSSIEN,ConvolutionR,length,width);
 }
 float diff(float a,float b){
 	

	if (a>b) return a-b;
	else return b-a;

}
/*Approximateur d'angles

params:
angle  : Angle a approximer à une valeur parmi (0°,45°,90°,135°)

return : 
ret    : valeur en rad de l'angle après approximation.

*/

float Approximer(float angle){
	


	
	if (angle >= 0.0  && angle < M_PI/4){
			
		if ( diff( angle , 0 ) < diff( angle , M_PI/4 ) ) { return 0.0;	}
		else return (M_PI/4);
	}

	else if (angle >= M_PI/4  && angle <= M_PI/2){
				
		if ( diff( angle , M_PI/4 ) < diff( angle , M_PI/2 ) ) { return M_PI/4;	}
		else return (M_PI/2);
	}

	


	else if (angle >= -M_PI/4  && angle < 0.0){
		if ( diff( angle , -M_PI/4 ) < diff( angle , 0.0 ) ) { return 3*M_PI/4;	}
		else return 0.0;
	}

	else if (angle >= -M_PI/2  && angle < -M_PI/4){
			
		if ( diff( angle , -M_PI/2 ) < diff( angle , -M_PI/4 ) ) { return M_PI/2;	}
		else return (3*M_PI/4);
	}

	
	return angle;
	}

/*fonction qui calcule pour ∀ (i,j) l'angle de la normale au gradient (approximé)
  et le module du gradient

  params: NormeGradient: matrice contenant la norme du gradient ∀ pt
  		  AngleGradient: matrice contenant les angles de la normale au gradient ∀pt
  		  ConvolutionR :      Image flouée
  		  length & width: mesures de toutes les matrices
*/
void NormeAngleGradient(float** NormeGradient,float ** AngleGradient, float** ConvolutionR, int length, int width ){  
    int dx,dy;
    float GradientEnX,GradientEnY;


  for(int i=0; i < length; i++){

    for(int j=0; j < width; j++){

      if ( j == 0 ) dx =  0;
      if ( i == 0 ) dy =  0;
      if ( j != 0 ) dx = j - 1 ;
      if ( i != 0 ) dy = i - 1 ;

      GradientEnX = ConvolutionR[i][dx] - ConvolutionR[i][j];
      GradientEnY = ConvolutionR[dy][j] - ConvolutionR[i][j];

      NormeGradient[i][j] = sqrt(GradientEnX * GradientEnX + GradientEnY * GradientEnY);


 	  AngleGradient[i][j] = (round)(Approximer(atan(-GradientEnX/ GradientEnY))*180/M_PI); 
 	   }
  }

  Recal(NormeGradient, length, width);
  SaveImagePgm(NAME_IMG_GRADIENT, NormeGradient, length, width);
}

/*function de suppression de non maximums
params: 
NormeGradient: Matrice contenant la norme du gradient ∀ pt
ResultSupp: NormeGradient apres la suppression des non maximums
*/

void SuppNonMaximums(float** ResultSupp,float** NormeGradient,float** AngleGradient,int length,int width){

	for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < width; j++)
			{
				ResultSupp[i][j]=NormeGradient[i][j];
			}
		}

		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < width; j++)
			{
				float ang=AngleGradient[i][j];


				if( ang == 0.0 && i - 1 > 0 && i + 1 < length - 1 && (NormeGradient[i][j] < NormeGradient[i-1][j] || NormeGradient[i][j] < NormeGradient[i+1][j]))
                	ResultSupp[i][j] = 0.0;


              	else if( ang == 45.0 && i - 1 > 0 && i + 1 < length - 1 &&  (NormeGradient[i][j] < NormeGradient[i][j-1] || NormeGradient[i][j] <  NormeGradient[i][j+1]))
              	{
                 	ResultSupp[i][j] = 0.0;
              	}
         
            	else if( ang == 90.0 && i - 1 > 0 && i + 1 < length - 1 && j - 1 > 0 && j + 1 < width - 1 &&( NormeGradient[i][j]  < NormeGradient[i - 1][j +1] || NormeGradient[i][j] < NormeGradient[i + 1][j - 1]))
                	ResultSupp[i][j] = 0.0;

              
               	else if( ang == 135.0 && i - 1 > 0 && i + 1 < length - 1 && j - 1 > 0 && j + 1 < width - 1 && (NormeGradient[i][j] <  NormeGradient[i + 1][j + 1] ||NormeGradient[i][j] <  NormeGradient[i - 1][j - 1]))
					ResultSupp[i][j] = 0.0;	

			}

		}

	Recal(ResultSupp, length, width);
	SaveImagePgm(NAME_IMG_SUPPRESSION, ResultSupp, length, width);
}

void follow(float** ResultHyst, float** ResultSupp, float** AngleGradient, float** dejacontour,int i, int j, int tau_L, int length, int width) {

    if(dejacontour[i][j]==1) return;//condition d'arret necessaire pour arreter la recursivité
        
   

   	dejacontour[i][j] = 1;	 	//le point est modifié ( pour eviter un appel infini de la recursivité)
    ResultHyst[i][j]  = 255.0; //mettre le point comme contour 

    
    if( AngleGradient[i][j]==0.0 ){
        if(i - 1 > 0 && ResultSupp[i - 1][j] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i - 1, j,tau_L, width, length);
        }
            
        if(i + 1 < length && ResultSupp[i + 1][j] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i + 1, j,tau_L, width, length);
        }
    }
        
 	if(AngleGradient[i][j]==45.0){
        if(i + 1 < length && j - 1 >= 0 && ResultSupp[i ][j - 1] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i , j - 1,tau_L, width, length);
        }
        
        if(j - 1 >= 0 && j + 1 < width && ResultSupp[i ][j + 1] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i , j + 1,tau_L, width, length);
        }
    }
       
        
    if(AngleGradient[i][j]==90.0){

        if(j - 1 >= 0 && j+1 >= width && i - 1 >= 0 && i+1 >= width  && ResultSupp[i-1 ][j + 1] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i-1 , j + 1,tau_L, width, length);
        }

        if(j - 1 >= 0 && j+1 >= width && i - 1 >= 0 && i+1 >= width && ResultSupp[i+1][j - 1] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i+1 , j - 1,tau_L, width, length);
        }
    	}

    if( AngleGradient[i][j] == 135.0){

        if(i - 1 >= 0 && j - 1 >= 0 && ResultSupp[i - 1][j - 1] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i - 1, j - 1, tau_L, width, length);
        }

        if(i + 1 < length && j + 1 < width && ResultSupp[i +  1][j + 1] >= tau_L) {
            follow(ResultHyst, ResultSupp, AngleGradient, dejacontour, i + 1, j + 1,tau_L, width,length);
        }
        
    	}
}
/*fonction de seuillage Hysteresis:
	pour tout point ayant une valeur > Th contour, on se dirige vers la direction
	de son gradient pour mettre les points ayant > Tl comme contour avec la fonction follow

  params: ResultHyst : matrice du resultat du seuillage par Hysteresis
  		  ResultSupp : matrice apres suppression des non-maximums
  		  AngleGradient: matrice contenant les angles de la normale au gradient ∀pt
  		  dejacontour: matrice [i][j]=1 le pixel est visité 0 sinon
  		  tau_h : borne maximale du seuil a ne pas depasser 
  		  tau_b : borne minimale du seuil a ne pas depasser 
  		  length & width: mesures de toutes les matrices
*/
void SeuillageHysteresis(float** ResultHyst, float** ResultSupp, float** AngleGradient, float** dejacontour, int tau_H,int tau_L, int length, int width){


	for(int  i = 0 ; i < length; i++ )
      	for(int j = 0 ; j < width; j++ ){

         	if( ResultSupp[i][j] >= tau_H )	follow(ResultHyst,ResultSupp, AngleGradient, dejacontour,i,j,tau_L, length, width);
            }

    Recal(ResultHyst, length, width);
	SaveImagePgm(NAME_IMG_CANNY, ResultHyst, length, width);

	
}

/*fonction calculant l'histogramme normalisé, puis l'histogramme cumulé (fct de répartition)
  on parcours l'histogramme cumulé et on recupere la plus grande valeur de la norme de gradient a laquelle ph<7

  param :  ResultCanny matrice de resultat de cannysemi auto
*/

void CannySemiAuto(float** ResultCanny, float** NormeGradient, float** ResultSupp, float **AngleGradient/*,float** dejacontour, */,float p_H, float* histo, int length, int width){
	

	int   inc=0;
	float tau_L=0;
	float tau_H1=0;
	float fonctionRepartition[256];

	
	compute_histo(NormeGradient,length,width,histo);
	float** dejacontour = fmatrix_allocate_2d(length,width);

	/*calcul fonction repartition*/
	fonctionRepartition[0]=histo[0];
	for (int i = 0; i <= 255; i++){
		fonctionRepartition[i] = fonctionRepartition[i - 1] + histo[i - 1];
	}


	for (int i = 0; i < 256; i++){

		inc=i;
		if (fonctionRepartition[i]>p_H ) break;
	}
	
	tau_H1=fonctionRepartition[inc-1]*100;
  	
  	tau_L = 0.5 * tau_H1;
  
  
   for(int  i = 0 ; i < length; i++ )
      	for(int j = 0 ; j < width; j++ ){

         	if( ResultSupp[i][j] >= tau_H1 )	follow(ResultCanny,ResultSupp, AngleGradient, dejacontour,i,j,tau_L, length, width);
            }

	Recal(ResultCanny, length, width);
	SaveImagePgm(NAME_IMG_CANNY_AUTO, ResultCanny, length, width);
}

int main(int argc,int** argv)
 {
  int i,j,k,l;
  int length,width;
  float tau_L;
  float tau_H;
  float p_H;
  float sigma;

  /*allocation matrice pour l'image photographe*/
  float** imgINR=LoadImagePgm(NAME_IMG_IN, &length, &width);
  float** imgINI=fmatrix_allocate_2d(length, width);

  /*allocation matrice pour le filtre gaussien*/
  float** IMRGaussian = fmatrix_allocate_2d(length, width);
  float** IMIGaussian = fmatrix_allocate_2d(length, width);

  /*allocation matrice contenant le resultat du (photographe * filtre gaussien)*/
  float** ConvolutionR = fmatrix_allocate_2d(length, width);
  float** ConvolutionI = fmatrix_allocate_2d(length, width);


  /*allocation des matrices norme gradient et angle gradient*/
  float** NormeGradient = fmatrix_allocate_2d(length, width);
  float** AngleGradient = fmatrix_allocate_2d(length, width);

  /*allocation matrice qui contient le resultat de la suppression */
  float** ResultSupp= fmatrix_allocate_2d(length, width);

  /*allocation matrice qui contient le resultat du seuillage par Hysteresis */
  float** ResultHyst= fmatrix_allocate_2d(length, width);

  /*allocation matrice pour stocker les cases deja parcouru*/
  float** dejacontour=fmatrix_allocate_2d(length, width);

  /*histogramme de dimension 256= tous les niveaux de gris possibles*/
  float histo[256];

  /*Matrice de resultat de canny semi auto*/
  float** ResultCanny=fmatrix_allocate_2d(length, width);


 

  for (int i = 0; i < length; i++)
  {
  	for (int j = 0; j < width; j++)
  	{	
  		AngleGradient[i][j]=0.0;
  		ResultHyst[i][j]=0.0;
  		dejacontour[i][j]=0.0;
  	}
  }
	
  /* Entrer des valeurs */
  printf("Entrez la valeur de tau_L: ");
  scanf("%f",&tau_L);
  printf("Entrez la valeur de tau_H: ");
  scanf("%f",&tau_H);
  printf("Entrez l'ecart type du filter Gaussien: ");
  scanf("%f",&sigma);
  
  /* Implementer detecteur de Canny */

  FiltreGaussien(ConvolutionR,ConvolutionI,imgINR,imgINI,IMRGaussian,IMIGaussian,sigma,length,width);
  NormeAngleGradient(NormeGradient, AngleGradient, ConvolutionR, length, width );

  /* suppression des non maximums */
  SuppNonMaximums(ResultSupp,NormeGradient,AngleGradient,length,width);
  

  /*seuillage par Hysteresis*/
  SeuillageHysteresis(ResultHyst,ResultSupp,AngleGradient,dejacontour,tau_H,tau_L,length,width);



  printf("Entrez la valeur de p_H: ");
  scanf("%f",&p_H);

   /*Canny semi auto*/
   CannySemiAuto(ResultCanny,NormeGradient,ResultSupp,AngleGradient/*,dejacontour*/,p_H,histo,length,width);




  /*liberation memoire*/
  free_fmatrix_2d(imgINR);
  free_fmatrix_2d(imgINI);
  free_fmatrix_2d(IMRGaussian);
  free_fmatrix_2d(IMIGaussian);
  free_fmatrix_2d(ConvolutionR);
  free_fmatrix_2d(ConvolutionI);
  free_fmatrix_2d(NormeGradient);
  free_fmatrix_2d(AngleGradient);
  free_fmatrix_2d(ResultSupp);
  free_fmatrix_2d(ResultHyst);
  free_fmatrix_2d(dejacontour);
  free_fmatrix_2d(ResultCanny);
 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}


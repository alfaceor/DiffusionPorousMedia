#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <time.h>

#define NPOSSIB 24 //  2x[ 6 (1stNeighr)] + 12 (2ndNeigh)

#define FREE 0
#define HARD 1
#define BORD 2

#define   MULTIPLIER  1664525
#define   INCREMENT   1013904223
#define   MODULUS     4294967295 //4294967296 


#define INDX(I,J,K) (I*UcellxSize*UcellySize + J*UcellySize + K)



using namespace std;
int main(const int argc, char** argv){
  time_t tstart, tend;
  tstart  = time(0);
  
  int   Hradius = -1;
  int   Bradius = -1;
  
  string strHradius = "";
  string strBradius = "";


  //--- Read options ---//
  int auxC;
  while (1)
    {
      static struct option long_options[] =
        {
          {"Hradius",  required_argument, 0, 'H'},
          {"Bradius",  required_argument, 0, 'B'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      auxC = getopt_long (argc, argv, "HB",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (auxC == -1)
        break;

      switch (auxC)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'H':
          printf ("option --Hradius %s \n", optarg);
          Hradius    = atoi(optarg);
          strHradius = optarg;
          break;

        case 'B':
          printf ("option --Bradius %s \n", optarg);
          Bradius    = atoi(optarg);
          strBradius = optarg;
          break;

        case '?':
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
        }
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }

  if ( Hradius < 0 ){
    printf("--Hradius is missing \n");
    exit (1);
  }

  if ( Bradius < 0 ){
    printf("--Bradius is missing \n");
    exit (1);
  }

  
  // General strings for files
  string flnPrefix = "rw3D_Geometry__Hrad_"+strHradius+"__Brad_"+strBradius;
  string strAux = "";
  stringstream sstrm;
  sstrm.str("");
  printf("#random walk 3D v 01 %s \n", flnPrefix.c_str());

//============================================================//

  const int Radius          = Hradius; //10;
  const int UcellxSize      = (2*(Radius)+1); // + 1;      
  const int UcellySize      = (2*(Radius)+1); // + 1;      
  const int UcellzSize      = (2*(Radius)+1); // + 1;      


  printf("# UcellxSize = %d, UcellySize = %d, UcellzSize = %d\n", UcellxSize, UcellySize, UcellzSize);
  
  // DEFINE THE NETWORK
  int currentPos = FREE;
  int nextPos    = FREE;
  int curr2next  = 10*currentPos + nextPos; // JUST TO GET A NOTATION USING MACROS
  
  int *net = new int[UcellxSize*UcellySize*UcellzSize]; // Three dimensional net
  for (int i =0; i<UcellxSize; i++)
    for (int j =0; j<UcellySize; j++)
      for (int k =0; k<UcellzSize; k++)
        net[INDX(i,j,k)] = 0;   


  //int r2  =  Radius*Radius;
  //int r2b = (Radius+border)*(Radius+border);
  int r2  = Hradius*Hradius; 
  int r2b = Bradius*Bradius; //  (Radius+border)*(Radius+border);
  int r2test = 0;
  int r2h = (Hradius-1)*(Hradius-1);  // border of the hard sphere.

  // Circle 
  int cx = Radius;
  int cy = Radius;
  int cz = Radius;

  for (int i= 0 ; i < UcellxSize; i++){
    for (int j=0 ; j < UcellySize; j++){
      for (int k=0 ; k < UcellzSize; k++){
        // net[i][j] = FREE;
        r2test = (i - cx)*(i - cx) + (j - cy)*(j - cy) + (k-cz)*(k-cz);
        if       ( r2test < r2h  )    net[INDX(i,j,k)] = -1;
        else  if ( r2test <= r2  )    net[INDX(i,j,k)] = HARD;
        else  if ( r2test  > r2b )    net[INDX(i,j,k)] = FREE;
        else                          net[INDX(i,j,k)] = BORD;
      }
    }
  }

  // Save network in a file
  sstrm.str(""); sstrm << flnPrefix << "__HARDWALL.csv"; strAux = sstrm.str();
  FILE *ofileNet_HARDWALL = fopen(strAux.c_str(), "w");
  sstrm.str(""); sstrm << flnPrefix << "__BORDER.csv"; strAux = sstrm.str();
  FILE *ofileNet_BORDER   = fopen(strAux.c_str(), "w");
  //fprintf(ofileNet, "#type i j k\n");
  fprintf(ofileNet_HARDWALL, "#i j k\n");
  fprintf(ofileNet_BORDER, "#i j k\n");

  int indx = 0; 
  int indy = 0; 
  int indz = 0; 
  for (int i=0; i<UcellxSize; i++){
    for (int j=0; j<UcellySize; j++){
      for (int k=0; k<UcellzSize; k++){
        // Get index in the cell unit
        indx = i % UcellxSize;
        if ( indx < 0 ) indx = indx + UcellxSize;
        indy = j % UcellySize;
        if ( indy < 0 ) indy = indy + UcellySize;
        indz = k % UcellzSize;
        if ( indz < 0 ) indz = indz + UcellzSize;

        // Check the type of the cell
        if (net[INDX(indx, indy, indz)] == HARD)
          fprintf(ofileNet_HARDWALL, "%d %d %d\n", i, j, k);
        else if (net[INDX(indx, indy, indz)] == BORD)
          fprintf(ofileNet_BORDER, "%d %d %d\n", i, j, k);
      }
    }
  }
  fclose(ofileNet_HARDWALL);
  fclose(ofileNet_BORDER);

  return EXIT_SUCCESS;
  //return 0;
}

